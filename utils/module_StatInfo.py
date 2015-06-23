from __future__ import division
import re,sys,os

class QcStat(object):
   def __init__(self,infile):
      self.infile = infile
      self.__init_statInfo()
   
   def __init_statInfo(self):
      self.raw_reads = 0
      self.raw_bases = 0
      self.cln_reads = 0
      self.cln_bases = 0
      self.ER_l_rate = 0.0
      self.ER_r_rate = 0.0
      self.Q20_l_rate= 0.0
      self.Q20_r_rate= 0.0
      self.Q30_l_rate= 0.0
      self.Q30_r_rate= 0.0
      self.GC_l_rate = 0.0
      self.GC_r_rate = 0.0
      self.N_remove  = 0
      self.Q_remove  = 0
      self.A_remove  = 0
      self.A_trimed  = 0
      
   def read_infile(self):
      if os.path.isfile( self.infile ):
         f_infile = open( self.infile,"r" )
         l_line = f_infile.readlines()
         f_1    = l_line[1].split('\t')
         self.raw_reads = int(f_1[0])
         self.raw_bases = int(f_1[1])
         self.cln_reads = int(f_1[2])
         self.cln_bases = int(f_1[3])
                  
         if len( f_1[5].split(';') ) > 1:
            self.ER_l_rate, self.ER_r_rate  = [ float(i) for i in f_1[5].split(';') ]
            self.Q20_l_rate,self.Q20_r_rate = [ float(i) for i in f_1[6].split(';') ]
            self.Q30_l_rate,self.Q30_r_rate = [ float(i) for i in f_1[7].split(';') ]
            self.GC_l_rate, self.GC_r_rate  = [ float(i) for i in f_1[8].split(';') ]
         else:
            self.ER_l_rate  = [ float(i) for i in f_1[5].split(';') ]
            self.Q20_l_rate = [ float(i) for i in f_1[6].split(';') ]
            self.Q30_l_rate = [ float(i) for i in f_1[7].split(';') ]
            self.GC_l_rate  = [ float(i) for i in f_1[8].split(';') ]
            
         self.N_remove  = l_line[2].split()[-1]
         self.Q_remove  = l_line[3].split()[-1]
         self.A_remove  = l_line[4].split()[2]
         self.A_trimed  = l_line[4].split()[-1]
         f_infile.close()

class TophatStat(dict):
   
   """
      Here I define a class for all tophat results.
      Because it is a dictionary class, ( TophatStat(dict) ), we can use  self['infile'] to define the input-file
   in the class. In the text-book, we used to use object class (  like TophatStat(object) ), in this condition, self['infile'] is written as self.infile.
   """
   def __init__(self,infile):
      """
         Initiate when:
            tophat_stat_info = TophatStat( in_file )
      """
      self['infile'] = infile
      self.__init_statInfo()
      self.__init_ItemPattern()

   """
      Public functions, which could be used out of the class
   """
   def read_infile(self):
      '''
         Read self['infile'] by line.
         1. strip command: remove \\n in the end of each lines.
         2. __determine_item, determine which item to get the words for Regular-Expression.
      '''
      f_infile = open( self['infile'],"r" )
      for line in f_infile:
         line = line.strip('\n')
         item,value = self.__determine_item(line)
         self.__load_itemValue( item,value )
      f_infile.close()
      
      
   def output_TophatStat(self):
      self['statInfo']['concordantPair'] = self['statInfo']['mappedPair'] - self['statInfo']['disconcordantPair']
      out_h = "Infile"
      out_1 = "%s" % (self['infile'])
      for item in ['totalRead','mappedRead','mappedPair','mappingRate','mappingConcordantRate','concordantPair','disconcordantPair']:
         value = self['statInfo'][item]
         out_h += "\t%s" % (item )
         out_1  = self.__out_itemValue( item,value, out_1 )
      print out_h
      print out_1
   
   
   """
      Private functions, which could be used ONLY by functions inside the class
   """   
   def __out_itemValue( self,item,value,  out_info ):
      item_type = self['itemType'][item]
      if    item_type == "str":
         out_info += "\t%s" % (value)
      elif  item_type == "int":
         out_info += "\t%d" % (value)
      return out_info
   
   def __init_statInfo( self):
      self['statInfo'] = {}
      self['statInfo']['totalRead']=0;
      self['statInfo']['mappedRead']=0;
      self['statInfo']['mappedPair']=0;
      self['statInfo']['mappingRate']=0;
      self['statInfo']['mappingConcordantRate']=0;
      self['statInfo']['concordantPair']=0;
      self['statInfo']['disconcordantPair']=0;
   
   def __init_ItemPattern(self):
      self['itemPattern'] = {                                                 \
         'totalRead'             :  "Input\s+:\s+(\d+)",                      \
         'mappedRead'            :  "Mapped\s+:\s+(\d+)",                     \
         'mappedPair'            :  "Aligned pairs:\s+(\d+)",                 \
         'mappingRate'           :  "^(.+?) overall read mapping rate",       \
         'mappingConcordantRate' :  "^(.+?) concordant pair alignment rate",  \
         'disconcordantPair'     :  "(\d+) (.+) are discordant alignments",   \
      }
      self['itemType'] = {                \
         'totalRead'             : "int", \
         'mappedRead'            : "int", \
         'mappedPair'            : "int", \
         'mappingRate'           : "str", \
         'mappingConcordantRate' : "str", \
         'concordantPair'        : "int", \
         'disconcordantPair'     : "int", \
      }
      
   def __determine_item(self,line):
      out_item = None
      out_value= None
      for item in self['itemPattern']:
         pat = self['itemPattern'][ item ]
         pattern = re.compile( pat )
         match   = pattern.search( line )
         if match:
            out_item = item
            out_value= match.group(1)
      return out_item, out_value

   def __load_itemValue( self, item,value ):
      if item != None:
         item_type = self['itemType'][item]
         if    item_type == "str":
            self['statInfo'][item]  =     value
         elif  item_type == "int":
            self['statInfo'][item] += int(value)
            
            
            
class SpikeIn(object):
   def __init__(self, file_HTS_SpikeIn,file_ERCC_info="/datc/huboqiang/cir_dyj_V2/Database/ercc.info.xls"):
      self.file_Spike_In_HTS = file_HTS_SpikeIn
      self.file_ERCC_info    = file_ERCC_info

      self.ERCC_info = { 'len':{}, 'mol':{} }
      self.ERCC_count= {}
      self.ERCC_total= 0
      self.RGC_count = {}
      self.__load_ERCC_info()
      
   def __load_ERCC_info(self):
      f_file = open( self.file_ERCC_info,"r")
      f_file.readline()
      for line in f_file:
         line = line.strip('\n')
         f    = line.split()
         ERCC_id = f[0]
         self.ERCC_info['len'][ERCC_id] = int(f[1])
         self.ERCC_info['mol'][ERCC_id] = float(f[2])
      f_file.close()
   
   def load_HTS_file(self):
      f_HTS = open( self.file_Spike_In_HTS )
      for line in f_HTS:
         line = line.strip('\n')
         f    = line.split()
         gene_id= f[0]
         if gene_id[0:4] == "ERCC":
            self.ERCC_count[ gene_id ] = int(f[1])
            self.ERCC_total           += int(f[1])
         else:
            self.RGC_count[  gene_id ] = int(f[1])
      f_HTS.close()
   
      
      
      