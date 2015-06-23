#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import subprocess
import cPickle as pickle
import time
import logging

import RNA_v2.utils.module_running_jobs as m_jobs
import RNA_v2.settings.projpath         as m_proj
import RNA_v2.settings.scripts          as m_scpt

logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
    stream=sys.stderr
)

class DataBaseInit(m_scpt.Scripts):
    
    def __init__(self, ref, sam_RNAinfo, is_debug=1):
        super(DataBaseInit,self).__init__()
        
        s_idx = ".".join(sam_RNAinfo.split("/")[-1].split(".")[:-1])
        self.define_scripts(s_idx)

        self.ref = ref
        self.is_debug = is_debug
    
    def check_files(self):
        self.define_files(self.ref)
        
        logging.info('Begin checking input files.')
        
        if not os.path.isdir("%s/%s" % (self.Database, self.ref)):
            os.mkdir("%s/%s" % (self.Database, self.ref))

        logging.info('Input database files were all put in %s/%s.'          %\
            (self.Database, self.ref) )
        
        '''
            1. Check if fasta file exists.
        '''
        if not os.path.isfile("%s.fa" % self.bt_index_base):
            logging.info('Input fasta %s not find. Now download from UCSC'  %\
                (self.bt_index_base) )
            time.sleep(0.5)
            self.__get_ref_fasta()
            logging.info('%s.fa generation done!' % (self.bt_index_base) )

        '''
            2. Check if fasta file were indexed.
        '''        
        fasta_idx = "%s.rev.1.bt2" % (self.bt_index_base)
        if not os.path.isfile(fasta_idx):
            logging.info('Fasta were not indexed.')
            time.sleep(0.5)
            logging.info('Now build index using bowtie2.')
            self.__get_ref_index()
            logging.info('Building index done!')

        '''
            3. Check if refGene file exists.
        '''
        if not os.path.isfile(self.genome_gtf):
            logging.info('Genome GTF file were not found.')
            time.sleep(0.5)
            logging.info('Now download refGene file from UCSC.')
            self.__get_refGene()
            logging.info('Generate refGene done!')
            
        '''
            4. Check if repeatMask file exists.
        '''
        if not os.path.isfile(self.rmsk_bed):
            logging.info('RepeatMask file were not found.')
            time.sleep(0.5)
            logging.info('Now download rmsk file from UCSC.')
            self.__get_rmsk()
            logging.info('Generate RepeatMask done!')

        '''
            5. Check if gtf transcriptome were built.
        '''
        fasta_idx = "%s.refGene.rev.1.bt2" % (self.bt_index_base)
        if not os.path.isfile(fasta_idx):
            logging.info('Transcriptome fasta file were not found.')
            time.sleep(0.5)
            logging.info('Now generating transcriptome fasta file.')
            self.__get_transcriptome()
            logging.info('Generate transcriptome fasta done!')
            
    
    def __get_ref_fasta(self):
        sh_file       = "%s/db01.DownloadRef.sh"      % (self.scripts)
        sh_work_file  = "%s/db01.DownloadRef_work.sh" % (self.scripts)
        
        l_sh_info = self.db_01_DownloadRef()
        l_sh_work = []
        l_sh_work.append("sh %s %s" % (sh_file,self.ref))
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)

    def __get_ref_index(self):
        sh_file       = "%s/db02.RefIndex.sh"      % (self.scripts)
        sh_work_file  = "%s/db02.RefIndex_work.sh" % (self.scripts)
        
        l_sh_info = self.db_02_BuildRefIndex()
        l_sh_work = []
        l_sh_work.append("sh %s %s" % (sh_file,self.ref))
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)

    def __get_refGene(self):
        sh_file       = "%s/db03.RefGene.sh"      % (self.scripts)
        sh_work_file  = "%s/db03.RefGene_work.sh" % (self.scripts)
        
        l_sh_info = self.db_03_RefGene()
        l_sh_work = []
        l_sh_work.append("sh %s %s" % (sh_file,self.ref))
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)


    def __get_rmsk(self):
        sh_file       = "%s/db04.rmsk.sh"      % (self.scripts)
        sh_work_file  = "%s/db04.rmsk_work.sh" % (self.scripts)
        
        l_sh_info = self.db_04_rmsk()
        l_sh_work = []
        l_sh_work.append("sh %s %s" % (sh_file,self.ref))
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)

    def __get_transcriptome(self):
        sh_file       = "%s/db05.transcriptome.sh"      % (self.scripts)
        sh_work_file  = "%s/db05.transcriptome_work.sh" % (self.scripts)
        
        l_sh_info = self.db_05_Transcriptome()
        l_sh_work = []
        l_sh_work.append("sh %s %s" % (sh_file,self.ref))
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)
    
    
