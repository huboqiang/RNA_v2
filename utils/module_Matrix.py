from __future__ import division
import re,sys,os
import module_StatInfo as Stat
import scipy
import numpy as np

class Matrix_info(object):
   def __init__(self,infile,inf_column=1,in_dtype="int"):
      self.infile     = infile
      self.inf_column = inf_column
      self.in_dtype   = in_dtype
      
   def load_mat(self):
      self.colname= []
      self.matrix = []
      self.rowname= []
      f_infile = open( self.infile,"r" )
      head = f_infile.readline()
      head = head.strip('\n')
      f_h  = head.split()
      self.rowname = f_h[ self.inf_column: ]
      for line in f_infile:
         line = line.strip('\n')
         f    = line.split()
         gene = "\t".join( f[0:self.inf_column] )
         self.colname.append( gene )
         self.matrix.append( np.array( f[self.inf_column:],dtype=self.in_dtype ) )
      f_infile.close()
      
      self.matrix = np.array( self.matrix,dtype=self.in_dtype )
