#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import subprocess
import cPickle as pickle
import time

import RNA_v2.utils.module_running_jobs as m_jobs
import RNA_v2.settings.projpath         as m_proj
import RNA_v2.settings.scripts          as m_scpt

def make_dir(l_args):
    """docstring for make_dir"""
    if len(l_args) == 1:
        """ only dictionary """
        if not os.path.isdir( l_args[0] ):
            os.mkdir( l_args[0] )
    
    elif len(l_args) == 2:
        """ dictionary/sample """
        if not os.path.isdir( l_args[0] ):
            os.mkdir( l_args[0] )
        if not os.path.isdir( "%s/%s" % (l_args[0],l_args[1]) ):
            os.mkdir(         "%s/%s" % (l_args[0],l_args[1]) )

    elif len(l_args) == 3:
        """ dictionary/sample/subname """
        if not os.path.isdir( l_args[0] ):
            os.mkdir( l_args[0] )
        if not os.path.isdir( "%s/%s"    % (l_args[0],l_args[1]) ):
            os.mkdir(         "%s/%s"    % (l_args[0],l_args[1]) )
        if not os.path.isdir( "%s/%s/%s" % (l_args[0],l_args[1],l_args[2]) ):
            os.mkdir(         "%s/%s/%s" % (l_args[0],l_args[1],l_args[2]) )


class Map_From_raw(m_scpt.Scripts):
    
    def __init__(self, ref, sam_RNAinfo, is_debug=1):
        super(Map_From_raw, self).__init__()
        
        self.s_idx = ".".join(sam_RNAinfo.split("/")[-1].split(".")[:-1])
        
        self.load_RNA_samInfo(sam_RNAinfo)
        self.define_scripts(self.s_idx)
        
        self.ref = ref
        self.is_debug = is_debug
        self.define_files(self.ref)
        self.M_CvtEnd = {"PE":2,"SE":1}
        
    def s01_QC(self):
        sh_file      = "%s/s01.QC.sh"      % (self.scripts)
        sh_work_file = "%s/s01.QC_work.sh" % (self.scripts)
        
        l_sh_info = self.s_01_QC()
        l_sh_work = []
        for samp in self.samInfo_pd_RNA['sample']:
            make_dir(  [ self.dir_clean_data, samp ] )
            idx       = (self.samInfo_pd_RNA['sample'] == samp)
            end       = self.samInfo_pd_RNA[ idx ]['end_type'].values[0]
            data_dype = self.M_CvtEnd[ end ]
            l_sh_work.append("sh %s %s %d" % (sh_file, samp, data_dype) )
      
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)


    def s02_Tophat(self):
        sh_file      = "%s/s02.Tophat.sh"      % (self.scripts)
        sh_work_file = "%s/s02.Tophat_work.sh" % (self.scripts)
        
        l_sh_info = self.s_02_Tophat()
        l_sh_work = []
        for samp in self.samInfo_pd_RNA['sample']:
            idx       =(self.samInfo_pd_RNA['sample'] == samp)
            brief_name= self.samInfo_pd_RNA[ idx ]['brief_name'].values[0]
            make_dir( [ self.dir_tophat, brief_name ] )
            end       = self.samInfo_pd_RNA[ idx ]['end_type'].values[0]
            data_dype = self.M_CvtEnd[ end ]
            l_sh_work.append("sh %s  %s %s %d"                              %\
                (sh_file, samp, brief_name, data_dype))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)

        