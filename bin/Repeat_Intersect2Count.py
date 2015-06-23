from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np

class Intersect(object):
   def __init__(self,infile):
      self.infile = infile
      self.pos_count = 0
      self.pos_infor = ""
      
   def load_file(self):
      f_infile = open( self.infile,"r" )
      
      self.__clean()
      for line in f_infile:
         line = line.strip('\n')
         f    = line.split()
         sub_group = f[9]
         group     = f[10].split("/")[0]
         element   = f[10].split("/")[-1]
         
         pos_info  = "%s\t%s\t%s\t%s" % ( "\t".join(f[0:3]), element, sub_group, group )
         
         if pos_info == self.pos_infor:
            self.__add( f[15] )
         else:
            if self.pos_infor != "":
               self.__print()
               self.__clean()
            self.__add( f[15] )
         self.pos_infor = pos_info
      
      if self.pos_infor != "":
         self.__print()
      
      f_infile.close()
         
         
   def __clean(self):
      self.pos_count = 0
      self.pos_infor = ""
      
   def __add(self,in_chrom):
      if in_chrom != ".":
         self.pos_count += 1
         
   def __print(self):
      print "%s\t%d" % ( self.pos_infor,self.pos_count )


def show_help():
   print >>sys.stderr,"\n\tpython",sys.argv[0],"/datc/huboqiang/tmp_lilin/PGC/repeat/Sample_PR10_TFC_141105-H9W-1.result2 "
   
def main():
   try:
      infile = sys.argv[1]
   except IndexError:
      show_help()
      sys.exit(1)

   samp_info = Intersect( infile )

   samp_info.load_file()






if __name__ == '__main__':
   main()
