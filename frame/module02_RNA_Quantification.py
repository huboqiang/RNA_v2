#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import subprocess
import cPickle as pickle
import time

import RNA_v2.utils.module_running_jobs as m_jobs
import RNA_v2.utils.module_CountInfo    as m_cnt
import RNA_v2.utils.module_GTFFeature   as m_gtf
import RNA_v2.settings.projpath         as m_proj
import RNA_v2.settings.scripts          as m_scpt
import RNA_v2.frame.module01_mapping_from_raw as m_01

class RNA_Quantification(m_scpt.Scripts):

    def __init__(self, ref, sam_RNAinfo, core_num=4, is_debug=1):
        super(RNA_Quantification, self).__init__()
        self.load_RNA_samInfo(sam_RNAinfo)
        self.s_idx = ".".join(sam_RNAinfo.split("/")[-1].split(".")[:-1])
        self.define_scripts(self.s_idx)

        self.ref = ref
        self.is_debug = is_debug
        self.core_num = 4
        self.define_files(self.ref)

    def run_pipeline(self, extra_GTF=None, given_GTF=None, is_MergeSam=0):
        """
                When not given a GTF file, this pipeline will try to detect
            noval transcripts. If a extra GTF file were given (extra means
            not included in UCSC refGene), it will be merged into a new gtf
            for gene quantification.
        """
        self.s03_HTSeq_known()
        if not given_GTF:
            self.s04_novo_1_cufflinks_u()
            self.s04_novo_2_cuffcomp_trans()
            self.s04_novo_3_HTSeq_unknown()
            self.s05_1_extra_makeGTF(extra_GTF)
            m_gtf.GTFFeature(self.merge_GTF)
            m_gtf.GTFFeature(self.merge_GTF_ercc)
            self.s06_1_cufflinks(self.dir_cufflinks_known, self.merge_GTF)
            self.s06_1_cufflinks(self.dir_cufflinks_known_ERCC, self.merge_GTF_ercc)
#            self.s06_cuffquant(self.dir_cuffquant, self.merge_GTF)
#            self.s06_cuffquant(self.dir_cuffquant_ERCC, self.merge_GTF_ercc)
            self.s08_repeatCount()
            if is_MergeSam:
                self.s07_1_mergeFPKM(
                    self.dir_cufflinks_known,
                    self.dir_cufflinks_known_mrg
                )
                self.s07_1_mergeFPKM(
                    self.dir_cufflinks_known_ERCC,
                    self.dir_cufflinks_known_ERCC_mrg
                )
#                self.s07_cuffnorm(
#                    self.dir_cuffquant_ERCC,
#                    self.dir_cuffnorm_ERCC,
#                    self.merge_GTF_ercc
#                )
#                self.s07_cuffnorm(
#                    self.dir_cuffquant,
#                    self.dir_cuffnorm,
#                    self.merge_GTF
#                )


        else:
            prefix = ".".join(given_GTF.split(".")[:-1])
            given_GTF_ERCC = "%s.ERCC.gtf" % (prefix)

            self.s05_2_extra_makeGTF(given_GTF, given_GTF_ERCC)
            m_gtf.GTFFeature(given_GTF)
            m_gtf.GTFFeature(given_GTF_ERCC)

            self.s06_1_cufflinks(self.dir_cufflinks_known, given_GTF)
            self.s06_1_cufflinks(self.dir_cufflinks_known_ERCC, given_GTF_ERCC)
#            self.s06_cuffquant(self.dir_cuffquant, given_GTF)
#            self.s06_cuffquant(self.dir_cuffquant_ERCC, given_GTF_ERCC)
            self.s08_repeatCount()
            if is_MergeSam:
                self.s07_1_mergeFPKM(
                    self.dir_cufflinks_known,
                    self.dir_cufflinks_known_mrg
                )
                self.s07_1_mergeFPKM(
                    self.dir_cufflinks_known_ERCC,
                    self.dir_cufflinks_known_ERCC_mrg
                )
#                self.s07_cuffnorm(
#                    self.dir_cuffquant,
#                    self.dir_cuffnorm,
#                    given_GTF
#                )
#                self.s07_cuffnorm(
#                    self.dir_cuffquant_ERCC,
#                    self.dir_cuffnorm_ERCC,
#                    given_GTF_ERCC
#                )

    def s03_HTSeq_known(self):
        sh_file      = "%s/s03.HTSeq.sh"      % (self.scripts)
        sh_work_file = "%s/s03.HTSeq_work.sh" % (self.scripts)

        l_sh_info = self.s_03_HTSeq_known()
        l_sh_work = []
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            m_01.make_dir([ self.dir_HTS_known, brief_name ])
            l_sh_work.append("sh %s  %s" % (sh_file, brief_name))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=self.core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)

    def s04_novo_1_cufflinks_u(self):
        sh_file      = "%s/s04.novo_1.cufflinks_unknown.sh"                 %\
            (self.scripts)
        sh_work_file = "%s/s04.novo_1.cufflinks_unknown_work.sh"            %\
            (self.scripts)

        l_sh_info = self.s_04_1_novo_HTSeq_unknown()
        l_sh_work = []
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            m_01.make_dir([ self.dir_cufflinks_unknown, brief_name ])
            l_sh_work.append("sh %s  %s" % (sh_file, brief_name))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=self.core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)


    def s04_novo_2_cuffcomp_trans(self):
        sh_file      = "%s/s04.novo_2.cuffcomp_unknown.sh"      % self.scripts
        sh_work_file = "%s/s04.novo_2.cuffcomp_unknown_work.sh" % self.scripts

        l_sh_info = self.s_04_2_novo_cuffcomp()
        l_sh_work = []

        dir_db = self.home_dir #"%s/%s" % (self.Database, self.ref)
        out_prefix = "%s/novo_lnc_raw_%s" % (dir_db, self.s_idx)
        l_input = []
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            infile    = "%s/%s/transcripts.gtf"                             %\
                (self.dir_cufflinks_unknown, brief_name)
            l_input.append(infile)

        insam = " ".join(l_input)
        l_sh_work.append("sh %s  %s %s" % (sh_file, out_prefix, insam))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=self.core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)

    def s04_novo_3_HTSeq_unknown(self):
        sh_file      = "%s/s04.novo_3.HTSeq_unknown.sh"      % self.scripts
        sh_work_file = "%s/s04.novo_3.HTSeq_unknown_work.sh" % self.scripts

        l_sh_info = self.s_04_3_novo_HTSeq_unknown()
        l_sh_work = []

        dir_db = self.home_dir #"%s/%s" % (self.Database, self.ref)
        unknown_GTF = "%s/novo_lnc_raw_%s.combined.gtf" % (dir_db, self.s_idx)
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            m_01.make_dir([ self.dir_HTS_unknown, brief_name ])
            l_sh_work.append("sh %s %s %s" %(sh_file, brief_name,unknown_GTF))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=self.core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)

        m_01.make_dir([ self.dir_HTS_result ])
        Gtf_Info = m_gtf.GTFFeature( unknown_GTF  )
        Gtf_Info.gene_intergenic( self.intragenic_bed )
        Cnt_Info = m_cnt.CountInfo(
            self.dir_HTS_unknown,
            self.samInfo_pd_RNA['brief_name'],
            "dexseq_NeoRaw",
            self.dir_HTS_result
        )

        Cnt_Info.generate_mat()
        Cnt_Info.load_mat()
        Cnt_Info.cal_RPKM(Gtf_Info.gene, self.dir_tophat)

        rpkm_file = "%s/merge.%s.RPKM.xls" % ( self.dir_HTS_result,
            "dexseq_NeoRaw"
        )

        Gtf_Info.load_gene_RPKM( rpkm_file )
        Gtf_Info.output_GTF()
        Gtf_Info.get_gene_info()

    def s05_1_extra_makeGTF(self, extra_GTF=None):
        sh_file      = "%s/s05.1.mergeGTF.sh"      % (self.scripts)
        sh_work_file = "%s/s05.1.mergeGTF_work.sh" % (self.scripts)

        l_sh_info = self.s_05_1_mergeGTF()
        l_sh_work = []

        dir_db = self.home_dir #"%s/%s" % (self.Database, self.ref)
        tail_str = "combined.FPKM0.5_rep0.25.multiExon.gtf"
        unknown_GTF = "%s/novo_lnc_raw_%s.%s" % (dir_db, self.s_idx, tail_str)

        l_sh_work.append("sh %s  %s %s" % (sh_file, unknown_GTF, extra_GTF))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=self.core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)

    def s05_2_extra_makeGTF(self, given_GTF, given_GTF_ERCC):
        sh_file      = "%s/s05.2.GTF_ercc.sh"      % (self.scripts)
        sh_work_file = "%s/s05.2.GTF_ercc_work.sh" % (self.scripts)

        l_sh_info = self.s_05_2_erccGTF(given_GTF, given_GTF_ERCC)
        l_sh_work = []

        l_sh_work.append("sh %s  %s" % (sh_file, self.ref))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=self.core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)


    def s06_cuffquant(self, cflk_dir, use_gtf):
        sh_file      = "%s/s06.cfqt.sh"      % (self.scripts)
        sh_work_file = "%s/s06.cfqt_work.sh" % (self.scripts)

        l_sh_info = self.s_06_cfqt(cflk_dir)
        l_sh_work = []

        for brief_name in self.samInfo_pd_RNA['brief_name']:
            m_01.make_dir([ cflk_dir, brief_name ])
            l_sh_work.append("sh %s  %s %s" % (sh_file, brief_name, use_gtf))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=self.core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)


    def s07_cuffnorm(self, cflk_dir, cfnm_dir, use_gtf):
        sh_file      = "%s/s07.cfnm.sh"      % (self.scripts)
        sh_work_file = "%s/s07.cfnm_work.sh" % (self.scripts)

        l_sh_info = self.s_07_cfnm(cfnm_dir)
        l_sh_work = []

        m_01.make_dir([ cfnm_dir ])
        l_brief = self.samInfo_pd_RNA['brief_name']
        l_cxb = [ "%s/%s/abundances.cxb" % (cflk_dir, sam) for sam in l_brief]

        list_sam = ",".join(l_brief)
        list_cxb = " ".join(l_cxb)

        l_sh_work.append(
            "sh %s  %s %s %s" % (sh_file, list_sam, use_gtf, list_cxb)
        )

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=self.core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)

    def s06_1_cufflinks(self, cflk_dir, use_gtf):
        sh_file      = "%s/s06.1.cflk.sh"      % (self.scripts)
        sh_work_file = "%s/s06.1.cflk_work.sh" % (self.scripts)

        l_sh_info = self.s_06_1_cflk(cflk_dir)
        l_sh_work = []

        for brief_name in self.samInfo_pd_RNA['brief_name']:
            m_01.make_dir([ cflk_dir, brief_name ])
            l_sh_work.append("sh %s  %s %s" % (sh_file, brief_name, use_gtf))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=self.core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)

    def s07_1_mergeFPKM(self, cflk_dir, mergeFPKM_dir):
        FPKM_Info = m_cnt.CountInfo(
            cflk_dir,
            self.samInfo_pd_RNA['brief_name'],
            "FPKM",
            mergeFPKM_dir
        )
        FPKM_Info.generate_mat()


    def s08_repeatCount(self):
        sh_file      = "%s/s08.RepCnt.sh"      % (self.scripts)
        sh_work_file = "%s/s08.RepCnt_work.sh" % (self.scripts)

        l_sh_info = self.s_08_repeatCount()
        l_sh_work = []
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            m_01.make_dir([ self.dir_repeat_counts, brief_name ])
            l_sh_work.append("sh %s  %s" % (sh_file, brief_name))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=self.core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)
