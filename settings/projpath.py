#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import subprocess
import numpy   as np
import cPickle as pickle
import pandas  as pd

class LoadSamp(object):
    def __init__(self):
        super( LoadSamp,self ).__init__()
        
    def load_RNA_samInfo(self,RNA_infile):
        self.samInfo_RNA_infile  = RNA_infile
        self.samInfo_pd_RNA  = pd.read_csv(self.samInfo_RNA_infile , sep="\t")
        
      
class DirSystem(object):
    def __init__(self):
        super( DirSystem,self ).__init__()
        self.home_dir               = os.path.abspath('./')
        path_file                   = os.path.realpath(__file__)
        self.dir_raw_data           = "%s/00.0.raw_data"               % (self.home_dir)
        self.dir_clean_data         = "%s/00.1.clean_data"             % (self.home_dir)
        self.dir_tophat             = "%s/01.Tophat"                   % (self.home_dir)
        self.dir_HTS_result         = "%s/02.HTSeq_result"             % (self.home_dir)
        self.dir_HTS_known          = "%s/02.1.HTSeq_known"            % (self.home_dir)
        self.dir_HTS_unknown        = "%s/02.2.HTSeq_unknown"          % (self.home_dir)
        self.dir_cufflinks_unknown  = "%s/03.cufflinks_unknown"        % (self.home_dir)
        self.dir_cuffquant          = "%s/04.cuffquant"                % (self.home_dir)
        self.dir_cuffnorm           = "%s/05.cuffnorm"                 % (self.home_dir)
        self.dir_cuffquant_ERCC     = "%s/04.1.cuffquant.ERCC"         % (self.home_dir)
        self.dir_cuffnorm_ERCC      = "%s/05.1.cuffnorm.ERCC"          % (self.home_dir)
        
        self.dir_cufflinks_known       = "%s/04.2.cufflinks_known"     % (self.home_dir)
        self.dir_cufflinks_known_ERCC  = "%s/04.2.cufflinks_known.ERCC"% (self.home_dir)
        
        self.dir_cufflinks_known_mrg       = "%s/05.2.cufflinksMerge_known"     % (self.home_dir)
        self.dir_cufflinks_known_ERCC_mrg  = "%s/05.2.cufflinksMerge_known.ERCC"% (self.home_dir)

        
        self.dir_repeat_counts      = "%s/06.repeat_counts"            % (self.home_dir)
        self.dir_repeat_mrg         = "%s/07.merge_repFPKM"            % (self.home_dir)
        self.dir_StatInfo           = "%s/StatInfo"                    % (self.home_dir)
        self.path                   = "%s"     % ("/".join(path_file.split('/')[:-2]))
        self.bin                    = "%s/bin" % ("/".join(path_file.split('/')[:-2]))

        if not os.path.exists( self.dir_StatInfo ):
            os.mkdir( self.dir_StatInfo )
        
        #####################################
        ######## Revise this path!!! ########
        #####################################
        self.Database               = "/data/Analysis/huboqiang/Database_RNA_v2"

        
    def define_scripts(self, s_idx):
        dir_script = "%s/scripts"      % (self.home_dir)
        self.scripts = "%s/scripts/%s" % (self.home_dir, s_idx)

        if not os.path.exists(dir_script):
            os.mkdir(dir_script)

        if not os.path.exists(self.scripts):
            os.mkdir(self.scripts)

        

class UsedSoftware(object):
    def __init__(self):
        super( UsedSoftware,self ).__init__()
        
        ##############################################
        ######## Revise the following path!!! ########
        ##############################################
        self.sftw_py        = "/home/huboqiang/anaconda/bin/python"
        self.sftw_pl        = "/data/Analysis/huboqiang/lib/local_perl/bin/perl"
        self.sftw_tophat_dir= "/data/Analysis/huboqiang/software/tophat-2.0.12.Linux_x86_64"
        self.sftw_cflk_dir  = "/data/Analysis/huboqiang/software/cufflinks-2.2.1.Linux_x86_64"
        self.sftw_bowtie_dir= "/data/Analysis/huboqiang/software/Bowtie/bowtie2-2.1.0"
        self.sftw_ucsc_dir  = "/data/Analysis/huboqiang/software/UCSC"
        self.sftw_samtools  = "/data/Analysis/huboqiang/software/samtools-0.1.18/samtools"
        self.sftw_bedtools  = "/data/Analysis/huboqiang/software/bedtools-2.17.0/bin/bedtools"
        self.sftw_deseq     = "/home/huboqiang/anaconda/lib/python2.7/site-packages/HTSeq/scripts/count.py"
        self.sftw_bgzip     = "/usr/local/bin/bgzip"
        self.sftw_tabix     = "/usr/local/bin/tabix"


class ProjInfo(LoadSamp,DirSystem,UsedSoftware):
    def __init__(self):
        super( ProjInfo,self ).__init__()

