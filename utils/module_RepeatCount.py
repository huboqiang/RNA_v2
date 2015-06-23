from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np
import matplotlib
import subprocess,time
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import RNA_v2.utils.module_StatInfo as m_Stat
import RNA_v2.utils.module_Matrix   as m_mat

class RepeatCount(object):
    def __init__(self,root_dir,l_samp,out_dir):
        self.l_sample = l_samp
        self.l_file   = [ "%s/%s/repeat_count.bed" % ( root_dir,samp ) for samp in l_samp ]
        self.out_dir  = out_dir

    def generate_mat(self):
        if not os.path.isdir( self.out_dir ):
            os.mkdir( self.out_dir )
        out_file_gene    = "%s/merge.Repeat.Count.xls" % ( self.out_dir )
        f_out_file_gene  = open( out_file_gene,"w" )

        out_info    = "Chrom\tBeg\tEnd\tElement\tSubGroup\tGroup\t%s" % ( "\t".join( self.l_sample ) )
        print >>f_out_file_gene, out_info

        shell_info = " paste %s " % (" ".join(self.l_file))
        p=subprocess.Popen(shell_info,stdout=subprocess.PIPE,shell=True)
        for line in p.stdout:
            line = line.strip('\n')
            f    = line.split()
            info = "\t".join(f[0:6])
            for i,val in enumerate(f):
                if i % 7 == 6:
                    info += "\t%d" % (int(val))
            
            print >>f_out_file_gene, info
        
        p.stdout.close()  
        f_out_file_gene.close()
      
    def element_group_subgroup_FPKM_sum(self,np_align_reads):
        infile = "%s/merge.Repeat.Count.xls" % ( self.out_dir )
        
        m_matrix = m_mat.Matrix_info( infile, 6, "float" )
        m_matrix.load_mat()
        
        outfile_element  = "%s/merge.Repeat.SumCount.element.xls"  % ( self.out_dir )
        outfile_subgroup = "%s/merge.Repeat.SumCount.subgroup.xls" % ( self.out_dir )
        outfile_group    = "%s/merge.Repeat.SumCount.group.xls"    % ( self.out_dir )
        
        f_outfile_element  = open( outfile_element ,"w" )
        f_outfile_subgroup = open( outfile_subgroup,"w" )
        f_outfile_group    = open( outfile_group   ,"w" )
        
        print >>f_outfile_element ,"Element\t%s"    % ( "\t".join( m_matrix.rowname ) )
        print >>f_outfile_subgroup,"SubGroup\t%s"   % ( "\t".join( m_matrix.rowname ) )
        print >>f_outfile_group   ,"Group\t%s"      % ( "\t".join( m_matrix.rowname ) )
        
        l_repInfo  = m_matrix.colname
        np_element  = np.array( [ "%s" % ( inf.split('\t')[3] ) for inf in l_repInfo ],dtype="string" )
        np_subgroup = np.array( [ "%s" % ( inf.split('\t')[4] ) for inf in l_repInfo ],dtype="string" )
        np_group    = np.array( [ "%s" % ( inf.split('\t')[5] ) for inf in l_repInfo ],dtype="string" )
        
        for elem in sorted( set(np_element) ):
            idx = (np_element ==elem)
            sub_mat      = m_matrix.matrix[ idx,: ]
            sum_FPKM     = np.sum( sub_mat,axis=0 )
            sum_FPKM    /= np_align_reads 
            sum_FPKM    *= 1000000 # normalize to each 100000
            str_sum_FPKM = np.array( sum_FPKM,dtype="string" )
            print >>f_outfile_element, "%s\t%s" % ( elem, "\t".join(str_sum_FPKM) )
        f_outfile_element.close()

        for subgroup in sorted( set(np_subgroup) ):
            idx = (np_subgroup ==subgroup)
            sub_mat      = m_matrix.matrix[ idx,: ]
            sum_FPKM     = np.sum( sub_mat,axis=0 )
            sum_FPKM    /= np_align_reads 
            sum_FPKM    *= 1000000 # normalize to each 100000
            str_sum_FPKM = np.array( sum_FPKM,dtype="string" )
            print >>f_outfile_subgroup, "%s\t%s" % ( subgroup, "\t".join(str_sum_FPKM) )
        f_outfile_subgroup.close()
        
        
        for group in sorted( set(np_group) ):
            idx = (np_group ==group)
            sub_mat      = m_matrix.matrix[ idx,: ]
            sum_FPKM     = np.sum( sub_mat,axis=0 )
            sum_FPKM    /= np_align_reads 
            sum_FPKM    *= 1000000 # normalize to each 100000
            str_sum_FPKM = np.array( sum_FPKM,dtype="string" )
            print >>f_outfile_group, "%s\t%s" % ( group, "\t".join(str_sum_FPKM) )
        f_outfile_group.close()