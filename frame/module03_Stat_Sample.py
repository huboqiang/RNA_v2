from __future__ import division
import re,sys,os

import scipy.stats
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

import RNA_v2.utils.module_StatInfo     as Stat
import RNA_v2.utils.module_CountInfo    as m_cnt
import RNA_v2.utils.module_GTFFeature   as m_gtf
import RNA_v2.utils.module_running_jobs as m_jobs
import RNA_v2.utils.module_RepeatCount  as m_repcnt
import RNA_v2.settings.scripts          as m_scpt

class SampStat(m_scpt.Scripts):
    
    def __init__(self, ref, sam_RNAinfo, given_GTF=None, is_debug=1):
        super(SampStat, self).__init__()
        self.load_RNA_samInfo(sam_RNAinfo)
        self.s_idx = ".".join(sam_RNAinfo.split("/")[-1].split(".")[:-1])
        
        self.define_scripts(self.s_idx)
        
        self.ref = ref
        self.is_debug = is_debug
        self.define_files(self.ref)

        self.__load_ERCC_info()
        
        dir_db = "%s/%s" % (self.Database, ref)
        self.use_gtf      = "%s/all.exon.sort.gtf"      % (dir_db)
        self.use_ercc_gtf = "%s/all.exon.sort.ERCC.gtf" % (dir_db)
        self.given_GTF = given_GTF
        if given_GTF:
            prefix = ".".join(given_GTF.split(".")[:-1])
            given_GTF_ERCC = "%s.ERCC.gtf" % (prefix)
            self.use_gtf      = given_GTF
            self.use_ercc_gtf = given_GTF_ERCC
        

    def Basic_Stat(self):
        """
            Stat for QC, Tophat mapping, ERCC RGC count
        """
        out_file   = "%s/01.BasicInfo_QC_map_SpikeIn.xls" % (self.dir_StatInfo)
        f_out_file = open( out_file,"w" )
        
        l_info = [
            "Rename", "Brief_name", "Raw_Reads", "Clean_Reads",
            "Pre_Map_Reads", "Aligned_Reads", "HTSseq_Known_Reads",
            "HTSeq_Refseq_Reads", "HTSeq_lncRNA_Reads", "HTSeq_Neo_Reads",
            "RFP_Reads", "GFP_Reads", "CRE_Reads", "ERCC_Reads",
            "RFP_polyA", "GFP_polyA", "CRE_polyA", "ERCC_Mols"
        ]
        
        out_info = "\t".join(l_info)

        print >>f_out_file, out_info

        l_brief = self.samInfo_pd_RNA['brief_name']

        ###   Load HTSeq reads
        HTS_info = m_cnt.CountInfo(
            self.dir_HTS_known,
            l_brief,
            "dexseq_clean",
            self.dir_HTS_result
        )
        HTS_info.load_mat(gen_col=1, generate = 1)
        HTS_info.sam_tot_reads()
        
        ###   Split refseq reads into RefSeq gene, lncGenes and novo genes.
        self.__get_HTS_clean_split()
        
        ###   Load HTSeq refseq reads
        Refseq_info = m_cnt.CountInfo(
            self.dir_HTS_known,
            l_brief,
            "dexseq_clean_refseq",
            self.dir_HTS_result
        )
        Refseq_info.load_mat(gen_col=1, generate = 0)
        Refseq_info.sam_tot_reads()
        
        ###   Load HTSeq lncRNA reads
        lncRNA_info = m_cnt.CountInfo(
            self.dir_HTS_known,
            l_brief,
            "dexseq_clean_lncRNA",
            self.dir_HTS_result
        )
        lncRNA_info.load_mat(gen_col=1, generate = 0)
        lncRNA_info.sam_tot_reads()

        ###   Load HTSeq novo-gene reads
        if not self.given_GTF:
            NeoPass_info = m_cnt.CountInfo(
                self.dir_HTS_unknown,
                l_brief,
                "dexseq_NeoPass",
                self.dir_HTS_result
            )
            NeoPass_info.load_mat(gen_col=1, generate = 0)
            NeoPass_info.sam_tot_reads()
        
        
        """
           Load other information
        """
        for samp in self.samInfo_pd_RNA['sample']:
            idx        = (self.samInfo_pd_RNA['sample'] == samp)
            brief_name =  self.samInfo_pd_RNA[idx]['brief_name'].values[0]
            rename     =  self.samInfo_pd_RNA[idx]['rename'].values[0]
            data_type  =  self.samInfo_pd_RNA[ idx ]['end_type'].values[0]
            QC_log       =  "%s/%s/log" % (self.dir_clean_data, samp)
            
            Tophat_log   =  "%s/%s/align_summary.txt"                       %\
                (self.dir_tophat, brief_name)
            
            HTSeq_SpikeIn = "%s/%s/%s.dexseq_ERCC_RGCPloyA.txt"             %\
                (self.dir_HTS_known, brief_name, brief_name)
            
            QcStat_info  = Stat.QcStat(QC_log)
            MapStat_info = Stat.TophatStat(Tophat_log)
            SpikeIn_info = Stat.SpikeIn(HTSeq_SpikeIn, self.ercc_info_file)
            
            
            QcStat_info.read_infile()
            MapStat_info.read_infile()
            SpikeIn_info.load_HTS_file()
            
            pre_map_read = MapStat_info['statInfo']['totalRead']
            aligned_read = MapStat_info['statInfo']['mappedRead']
            
            if data_type == "PE":         
                HTSseq_read  = self.__get_HTS_reads(HTS_info    ,samp) * 2
                Refseq_read  = self.__get_HTS_reads(Refseq_info ,samp) * 2
                lncRNA_read  = self.__get_HTS_reads(lncRNA_info ,samp) * 2
                NeoPass_read = 0
                if self.given_GTF:
                    NeoPass_read = self.__get_HTS_reads(NeoPass_info,samp) * 2
                read_RFP     = SpikeIn_info.RGC_count['RGC-mRFP'] * 2
                read_GFP     = SpikeIn_info.RGC_count['RGC-GFP' ] * 2
                read_CRE     = SpikeIn_info.RGC_count['RGC-CRE' ] * 2
                read_ERCC    = SpikeIn_info.ERCC_total            * 2
            else:
                HTSseq_read  = self.__get_HTS_reads(HTS_info    ,samp)
                Refseq_read  = self.__get_HTS_reads(Refseq_info ,samp)
                lncRNA_read  = self.__get_HTS_reads(lncRNA_read ,samp)
                if self.given_GTF:
                    NeoPass_read = self.__get_HTS_reads(NeoPass_info,samp)

                read_RFP     = SpikeIn_info.RGC_count['RGC-mRFP']
                read_GFP     = SpikeIn_info.RGC_count['RGC-GFP' ]
                read_CRE     = SpikeIn_info.RGC_count['RGC-CRE' ]
                read_ERCC    = SpikeIn_info.ERCC_total
            
            mol_RFP  = self.samInfo_pd_RNA[idx]['RFP_polyA'].values[0]
            mol_GFP  = self.samInfo_pd_RNA[idx]['GFP_polyA'].values[0]
            mol_CRE  = self.samInfo_pd_RNA[idx]['CRE_polyA'].values[0]
            mol_ERCC = self.samInfo_pd_RNA[idx]['ERCC_time'].values[0]   *\
                        6.023*10**10
            
            l_out = [
                rename, brief_name,
                QcStat_info.raw_reads, QcStat_info.cln_reads,
                pre_map_read, aligned_read, HTSseq_read,
                Refseq_read,  lncRNA_read,  NeoPass_read,
                read_RFP, read_GFP,read_CRE, read_ERCC,
                mol_RFP , mol_GFP ,mol_CRE, mol_ERCC
            ]
            l_out = [str(i) for i in l_out]
            print  >>f_out_file, "\t".join(l_out)
            
        f_out_file.close()



    
    def __load_ERCC_info(self):
        f_file = open(self.ercc_info_file, "r")
        
        self.ERCC_info = { 'len':{}, 'mol':{} }
        f_file.readline()
        for line in f_file:
            line = line.strip('\n')
            f    = line.split()
            ERCC_id = f[0]
            self.ERCC_info['len'][ERCC_id] = int(f[1])
            self.ERCC_info['mol'][ERCC_id] = float(f[2])
        f_file.close()
        
    def __get_HTS_reads(self, HTS_info, samp):
        val = 0
        if samp in self.samInfo_pd_RNA['sample']:
            idx = list(self.samInfo_pd_RNA['sample']).index( samp )
            val = HTS_info.sam_tot_reads[idx]
        return val
        
    def __get_HTS_clean_split(self):
        sh_file      = "%s/stat02.SplitHTS.sh"      % (self.scripts)
        sh_work_file = "%s/stat02.SplitHTS_work.sh" % (self.scripts)
        
        l_sh_info = self.stat_02_SplitHTS()
        l_sh_work = []
            
        l_sh_work.append("sh %s %s" % (sh_file, self.use_gtf))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)





    def ERCC_Stat(self):
        """
            Using ERCC mols to estimate the total amount of 
            refGene RNA-mols in a cell. 
            (ERCC_MOLs) / (ERCC_FPKM) = (mRNA_MOLs) / (mRNA_FPKM)
        """
        self.l_ERCC_name = []
        self.l_RGCs_name = []
        self.l_mRNA_name = []
        self.l_ERCC_FPKM = {}
        self.l_RGCs_FPKM = {}
        self.l_mRNA_FPKM = {}
        self.l_cirRNA_FPKM={}
        
        self.l_ERCC_HTSname = []
        self.l_RGCs_HTSname = []
        self.l_mRNA_HTSname = []
        self.l_ERCC_RPKM = {}
        self.l_RGCs_RPKM = {}
        self.l_mRNA_RPKM = {}
        
        
        self.l_ERCC_MOLs = {}
        self.l_RGCs_MOLs = {}
        self.l_mRNA_MOLs = {}
        self.l_cirRNA_MOLs={}
        self.l_mRNA_MOLs_HTSname = {}
        
        self.regression = {}
        
        self.__load_FPKM()
        self.__load_MOLs()      # ERCC RGC mols
        self.__get_mRNA_MOLs()  # get mRNA mols using ERCC_FPKM, ERCC_MOLs and mRNA_FPKM
        self.__load_Count()
        
        out_file   = "%s/02.ERCC_Mols.xls" % (self.dir_StatInfo)
        f_out_file = open( out_file,"w" )
        
        l_info = [
            "Sample", "Brief_samp", "ERCC_MOLs", "RGC_MOLs", "mRNA_MOLs",
            "RefSeq_mRNA_MOLs", "Regression_R",
            "Regression_P", "RefSeq_mRNA_FPKM>0.1"
        ] 
        print >>f_out_file, "\t".join(l_info)
        
        for samp in self.samInfo_pd_RNA['sample']:
            idx        = (self.samInfo_pd_RNA['sample'] == samp)
            brief_name =  self.samInfo_pd_RNA[idx]['brief_name'].values[0]
            rename     =  self.samInfo_pd_RNA[idx]['rename'].values[0]

            ERCC_MOLs  = sum(   self.l_ERCC_MOLs[brief_name])
            RGC_MOLs   = sum(   self.l_RGCs_MOLs[brief_name])
            mRNA_MOLs  = np.sum(self.l_mRNA_MOLs[brief_name])
            RefSeq_mRNA_MOLs  =                                              \
                np.sum(self.l_mRNA_MOLs[brief_name][self.mRNA_refSeq_index])
                
            RefSeq_mRNA_lFPKM =                                              \
                np.array(self.l_mRNA_FPKM[brief_name],dtype=float)
                
            RefSeq_mRNA_lFPKM = RefSeq_mRNA_lFPKM[ self.mRNA_refSeq_index ]
            RefSeq_mRNA_Exps  =                                              \
                np.shape(RefSeq_mRNA_lFPKM[RefSeq_mRNA_lFPKM >= 0.1])[0]
                
            regression_R      = self.regression[brief_name]['r_value']
            regression_P      = self.regression[brief_name]['p_value']
            
            l_out = [
                rename, brief_name,
                ERCC_MOLs, RGC_MOLs, mRNA_MOLs, RefSeq_mRNA_MOLs,
                regression_R, regression_P, RefSeq_mRNA_Exps
            ]
            l_out = [str(i) for i in l_out]
            print >>f_out_file, "\t".join(l_out)
        f_out_file.close()
        
    def __load_Count(self):  # for ERCC count Stat
        l_brief = self.samInfo_pd_RNA['brief_name']
        Gtf_Info = m_gtf.GTFFeature(self.use_ercc_gtf)
        ERCC_info = m_cnt.CountInfo(
            self.dir_HTS_known,
            l_brief,
            "dexseq_ERCC_RGCPloyA",
            self.dir_HTS_result
        )
        ERCC_info.load_mat(gen_col=1, generate=1)
        ERCC_info.cal_RPKM(Gtf_Info.gene,self.dir_tophat)
        
    def __load_FPKM(self):
        Cuffnorm_file = "%s/merge.FPKM.gene.xls" % ( self.dir_cufflinks_known_ERCC_mrg )
        f_Cuffnorm    = open( Cuffnorm_file,"r" )
        in_head       = f_Cuffnorm.readline()
        f_head        = in_head.split()
        l_samp_brief  = [ sam for sam in f_head[1:] ]
        
        for brief_name in l_samp_brief:
            self.l_ERCC_FPKM[brief_name] = []
            self.l_RGCs_FPKM[brief_name] = []
            self.l_mRNA_FPKM[brief_name] = []
        
        for line in f_Cuffnorm:
            line = line.strip('\n')
            f    = line.split()
            gene = f[0]
            
            if gene[0:3].upper() == "MIR":
                continue
            if gene[0:5].upper() == "SNORD":
                continue
            
            if gene[0:5].upper() == "ERCC-":
                self.l_ERCC_name.append(gene)
                for i,brief_name in enumerate(l_samp_brief):
                    self.l_ERCC_FPKM[brief_name].append(f[i+1])
                 
            elif gene[0:4].upper() == "RGC-":
                self.l_RGCs_name.append(gene)
                for i,brief_name in enumerate(l_samp_brief):
                    self.l_RGCs_FPKM[brief_name].append(f[i+1])
            
            else:
                self.l_mRNA_name.append(gene)
                for i,brief_name in enumerate(l_samp_brief):
                    self.l_mRNA_FPKM[brief_name].append(f[i+1])
        
        f_Cuffnorm.close()
        self.__mRNA_refSeq_index()
        
    def __mRNA_refSeq_index(self):
        self.mRNA_refSeq_index = []
        for i,gene in enumerate(self.l_mRNA_name):
            if  gene[0:7].upper() == "lncGene" or \
                gene[0:9].upper() == "novoXLOC_":
                continue
            self.mRNA_refSeq_index.append( i )
        self.mRNA_refSeq_index = np.array( self.mRNA_refSeq_index,dtype=int )
    
  
    def __load_MOLs(self):
        for samp in self.samInfo_pd_RNA['sample']:
            idx        = (self.samInfo_pd_RNA['sample'] == samp)
            brief_name =  self.samInfo_pd_RNA[idx]['brief_name'].values[0]
            rename     =  self.samInfo_pd_RNA[idx]['rename'].values[0]
            
            if brief_name not in self.l_ERCC_MOLs:
                self.l_ERCC_MOLs[brief_name] = []
            if brief_name not in self.l_RGCs_MOLs:
                self.l_RGCs_MOLs[brief_name] = []
               
            for ERCC_id in self.l_ERCC_name:

                ERCC_time = self.samInfo_pd_RNA[idx]['ERCC_time'].values[0]
                self.l_ERCC_MOLs[brief_name].append(
                    self.ERCC_info['mol'][ERCC_id]*(6.02*10**23/10**18)*ERCC_time
                )
            

            
            CRE_polyA = self.samInfo_pd_RNA[idx]['CRE_polyA'].values[0]
            GFP_polyA = self.samInfo_pd_RNA[idx]['GFP_polyA'].values[0]
            RFP_polyA = self.samInfo_pd_RNA[idx]['RFP_polyA'].values[0]
            
            self.l_RGCs_MOLs[brief_name].append(CRE_polyA)
            self.l_RGCs_MOLs[brief_name].append(GFP_polyA)
            self.l_RGCs_MOLs[brief_name].append(RFP_polyA)
         
    def __get_mRNA_MOLs(self):
        self.regression = {}
        for samp in self.samInfo_pd_RNA['sample']:
            idx        = (self.samInfo_pd_RNA['sample'] == samp)
            brief_name =  self.samInfo_pd_RNA[idx]['brief_name'].values[0]
            rename     =  self.samInfo_pd_RNA[idx]['rename'].values[0]
#            print self.l_ERCC_FPKM
#            print self.l_ERCC_FPKM.keys()
            np_ERCC_FPKM = np.array(self.l_ERCC_FPKM[brief_name],dtype=float)
            np_ERCC_MOLs = np.array(self.l_ERCC_MOLs[brief_name],dtype=float)
            
            np_ERCC_FPKM = np.log2( np_ERCC_FPKM+0.1 )
            np_ERCC_MOLs = np.log2( np_ERCC_MOLs )
            
            mol_idx = (np_ERCC_MOLs-np.log2(6.02*10**23/10**18) > -18) *\
                (np_ERCC_FPKM>0)
            
            np_ERCC_FPKM = np_ERCC_FPKM[mol_idx]
            np_ERCC_MOLs = np_ERCC_MOLs[mol_idx]
            
            slope,intercept,r_value,p_value,slope_std_error                 =\
                scipy.stats.linregress(np_ERCC_MOLs, np_ERCC_FPKM)
            
            self.regression[brief_name] = {
                'slope'    :slope,
                'inter'    :intercept,
                'r_value'  :r_value,
                'p_value'  :p_value,
                'std_err'  :slope_std_error
            }
            
            np_mRNA_FPKM = np.array(self.l_mRNA_FPKM[brief_name],dtype="float")
            np_mRNA_FPKM = np.log2(np_mRNA_FPKM + 0.1)
            np_mRNA_MOLs = (np_mRNA_FPKM - intercept) / slope
            np_mRNA_MOLs = np.power(np_mRNA_MOLs, 2)
            np_mRNA_MOLs[np_mRNA_MOLs<0.0]           = 0.0
            np_mRNA_MOLs[np_mRNA_FPKM<=np.log2(0.1)] = 0.0
            
            self.l_mRNA_MOLs[brief_name] = np_mRNA_MOLs
            


    def Repeat_Stat(self):
        if not os.path.isdir(self.dir_repeat_mrg):
            os.mkdir(self.dir_repeat_mrg)
            
        l_brief = self.samInfo_pd_RNA['brief_name']
        
        l_align_reads = []
        for brief_samp in l_brief:
            prefix       = "%s/%s" % (self.dir_tophat, brief_samp)
            Tophat_log   = "%s/align_summary.txt" % (prefix)
            MapStat_info = Stat.TophatStat(Tophat_log)
            MapStat_info.read_infile()
            l_align_reads.append(MapStat_info['statInfo']['mappedRead'])
           
        np_align_reads = np.array(l_align_reads)
        
        RepCnt = m_repcnt.RepeatCount(
            self.dir_repeat_counts, 
            l_brief, 
            self.dir_repeat_mrg
        )
        RepCnt.generate_mat()
        RepCnt.element_group_subgroup_FPKM_sum(np_align_reads)