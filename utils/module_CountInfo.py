from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np
import matplotlib
import subprocess,time

import RNA_v2.utils.module_StatInfo as m_Stat

class CountInfo(object):
    def __init__(self,root_dir,l_samp,ltype,out_dir):
        self.l_sample = l_samp
        self.l_file   = [ 
            "%s/%s/%s.%s.txt" % (root_dir,samp,samp,ltype) for samp in l_samp
        ]
        self.out_dir  = out_dir
        self.ltype    = ltype


    def generate_mat(self):
        if not os.path.isdir( self.out_dir ):
            os.mkdir( self.out_dir )
        out_file_gene    = "%s/merge.%s.gene.xls" % (self.out_dir, self.ltype)
        f_out_file_gene  = open(out_file_gene,"w" )
        out_file_stat    = "%s/merge.%s.stat.xls" % (self.out_dir, self.ltype)
        f_out_file_stat  = open(out_file_stat,"w")
        
        out_info    = "Gene\t%s" % ("\t".join( self.l_sample))
        print >>f_out_file_gene, out_info
        print >>f_out_file_stat, out_info
        
        shell_info = "paste %s" % (" ".join( self.l_file))     
        p=subprocess.Popen(shell_info,stdout=subprocess.PIPE,shell=True)
        for line in p.stdout:
            line = line.strip('\n')
            f    = line.split()
            info = f[0]
            for i,val in enumerate(f[1:]):
                if i % 2 == 0:
                    info += "\t%s" % (val)
            
            if line[0:2] != "__":
                print >>f_out_file_gene, info
            else:
                print >>f_out_file_stat, info
        
        p.stdout.close()  
        f_out_file_gene.close()
        f_out_file_stat.close()
    
    def load_mat(self,infile="",gen_col=1, generate = 0):
        self.gene = []
        self.cnt_matrix = []
        
        if infile == "":
            infile = "%s/merge.%s.gene.xls"%(self.out_dir,self.ltype)

        if generate:
            self.generate_mat()
        
        f_infile = open(infile,"r")
        f_infile.readline()
        for line in f_infile:
            line = line.strip('\n')
            f    = line.split()
            gene = ",".join(f[0:gen_col])
            self.gene.append(gene)
            self.cnt_matrix.append(np.array(f[gen_col:],dtype=int))
        f_infile.close()
        
        self.cnt_matrix = np.array(self.cnt_matrix,dtype=int)
    
    def sam_tot_reads(self):
        self.sam_tot_reads = np.sum(self.cnt_matrix,axis=0)
    
    def cal_RPKM(self,M_gen_len,tophat_dir):
        self.RPKM      = []
        self.RPKM_gene = []
        l_sam_TotReads = []
        for samp in self.l_sample:
            align_log   = "%s/%s/align_summary.txt" % (tophat_dir,samp)
            StatInfo    = m_Stat.TophatStat( align_log )
            StatInfo.read_infile()
            l_sam_TotReads.append(StatInfo['statInfo']['mappedPair'])
        
        np_sam_TotReads = np.array( l_sam_TotReads,dtype=int )
        for i,gene in enumerate(self.gene):
            if  (gene not in M_gen_len) or                                   \
                ('max_len' not in M_gen_len[gene]) or                        \
                (M_gen_len[gene] == 0):
                print >>sys.stderr,                                          \
                   "%s not in Gene_len list or with not valid maxlen" % (gene)
                continue
            self.RPKM_gene.append(gene)
            max_len = M_gen_len[gene]['max_len']
            l_RPKM  = (self.cnt_matrix[i]*10**9)/(max_len*np_sam_TotReads)
            self.RPKM.append(l_RPKM)
        self.RPKM = np.array(self.RPKM,dtype=float)
        
        self.__write_RPKM()
    
    
    def __write_RPKM(self):
        out_file_RPKM    = "%s/merge.%s.RPKM.xls" % (self.out_dir, self.ltype)
        f_out_file_RPKM  = open(out_file_RPKM,"w")
        out = "Gene\t%s" % ("\t".join( self.l_sample))
        print >>f_out_file_RPKM, out
        for i,gene in enumerate(self.RPKM_gene):
            lis_RPKM = "\t".join(np.array(self.RPKM[i],dtype="string"))
            print >>f_out_file_RPKM, "%s\t%s" % (gene,lis_RPKM)
        f_out_file_RPKM.close()
