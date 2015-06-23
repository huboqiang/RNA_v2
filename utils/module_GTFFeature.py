from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np
import matplotlib
import subprocess,time
matplotlib.use('Agg')
from matplotlib import pyplot as plt

class GTFFeature(object):
   def __init__(self, gtf_file):
      self.gtf       =  gtf_file
      self.tran   =  {}
      self.gene   =  {}
      self.__load_GTF()
      self.RPKMInfo = {}
      
   def __load_GTF(self):

      pat_gene = "gene_id\s+\"(\S+)\";"
      pat_tran = "transcript_id\s+\"(\w+)\";"
      pat_exon = "exon_number\s+\"*(\w+)\"*"

      pattern_gene = re.compile( pat_gene )
      pattern_tran = re.compile( pat_tran )
      pattern_exon = re.compile( pat_exon )

      f_infile = open( self.gtf, "r" )
      for line in f_infile:
         line  =  line.strip('\n')
         f     =  line.split('\t')
         
         match_gene = pattern_gene.search( f[-1] )
         match_tran = pattern_tran.search( f[-1] )
         match_exon = pattern_exon.search( f[-1] )
         
         gene = ""
         tran = ""
         exon = ""
         if match_gene:
            gene = match_gene.group(1)
         if match_tran:
            tran = match_tran.group(1)
         if match_exon:
            exon = match_exon.group(1)
         
         if gene == "" or tran == "" or exon == "":
            if gene[0:5] != "ERCC-" and gene[0:4] != "RGC-":
               print  "ERROR line: %s; Gene %s, Tran %s, Exon %s" % (line, gene,tran,exon )
            else:
               exon = "1"
         
         if gene not in self.gene:
            self.gene[ gene ] = { 'tran':[], 'max_len':0, 'exon':{}, 'line':{}, 'is_intergenic':0  }
            
         if tran not in self.gene[gene]['tran']:
            self.gene[gene]['tran'].append( tran )
            
         if tran not in self.gene[gene]['exon']:
            self.gene[gene]['exon'][tran] = []
         self.gene[gene]['exon'][tran].append( exon )
         
         if tran not in self.gene[gene]['line']:
            self.gene[gene]['line'][tran] = {}
         self.gene[gene]['line'][tran][exon] = line

         tmp_len =   int(f[4]) - int(f[3])
         if tran not in self.tran:
            self.tran[ tran ] = { 'gene':gene, 'tmp_len':0, 'exon_cnt':0, 'beg':int(f[3]), 'end':int(f[4]), 'chr':f[0] }
         self.tran[ tran ][ 'tmp_len'  ] += tmp_len
         self.tran[ tran ][ 'exon_cnt' ] += 1
         if int(f[3]) < self.tran[ tran ]['beg']:
            self.tran[ tran ]['beg'] = int(f[3])
         if int(f[4]) > self.tran[ tran ]['end']:
            self.tran[ tran ]['end'] = int(f[4])

      f_infile.close()
   
      out_prefix = "%s"     % ( ".".join( self.gtf.split(".")[:-1] ) )
      out_file_bed = "%s.gene.bed" % ( out_prefix )
      f_out_file_bed = open( out_file_bed,"w" )
      out_file_len = "%s.genlen"   % ( out_prefix )
      f_out_file_len = open( out_file_len,"w" )

      
      # print gene length
      for gene in sorted( self.gene ):
         max_len = 0
         for tran in self.gene[gene]['tran']:
            if max_len < self.tran[ tran ][ 'tmp_len' ]:
               max_len = self.tran[ tran ][ 'tmp_len' ]
      
         self.gene[ gene ]['max_len'] =  max_len
         
         for tran in sorted( self.gene[gene]['tran'] ):
            print >>f_out_file_len,"%s\t%s\t%d\t%d\t%d" % (  gene,tran, self.tran[ tran ]['tmp_len'],max_len, self.tran[ tran ][ 'exon_cnt' ]  )
            
      # print transcript bed
         for tran in sorted( self.gene[gene]['tran'] ):
            print >>f_out_file_bed,"%s\t%d\t%d\t%s\t%s" % ( self.tran[tran]['chr'], self.tran[tran]['beg'], self.tran[tran]['end'],tran,gene )
      
      f_out_file_bed.close()
      f_out_file_len.close()
      
      
   def gene_intergenic( self,intragenic_bed ):
      
      out_prefix = "%s"     % ( ".".join( self.gtf.split(".")[:-1] ) )
      out_file_bed = "%s.gene.bed" % ( out_prefix )
      
      shell_info = "bedtools closest -d -a %s -b %s " % ( out_file_bed, intragenic_bed )
      p=subprocess.Popen(shell_info,stdout=subprocess.PIPE,shell=True)
      for line in p.stdout:
         line = line.strip('\n')
         f    = line.split()
         gene = f[4]
         dist = int(f[-1])
         if dist < 10000:
            continue
         else:
            self.gene[ gene ]['is_intergenic'] = 1
      p.stdout.close()
      
   def load_gene_RPKM(self,rpkm_file):
      self.RPKMInfo['samp'] = {}
      self.RPKMInfo['gene_info'] = {}
      self.RPKMInfo['list_gene'] = []
      
      print "opening %s" % (rpkm_file)
      f_sam_file_rpkm = open( rpkm_file,"r" )
      in_h = f_sam_file_rpkm.readline()
      self.RPKMInfo['samp']['list'] = in_h.split()[1:]
      
      self.__load_tag()
      for line in f_sam_file_rpkm:
         line = line.strip('\n')
         f    = line.split()
         gene = f[0]
         result = self.__is_pass( line )
         if gene not in self.RPKMInfo['gene_info']:
            self.RPKMInfo['gene_info'][gene] = { 'is_pass':result,'trans':[], 'exons':{},'line':{}  }
         self.RPKMInfo['gene_info'][gene]['is_pass'] = result
         self.RPKMInfo['list_gene'].append( gene )
      f_sam_file_rpkm.close()
   
   def output_GTF(self):
      self.use_gene = {}
      self.use_tids = {}
      
      gtf_file = self.gtf
      out_prefix = "%s" % ( ".".join( gtf_file.split(".")[:-1] ) )
      f_out = open( "%s.FPKM0.5_rep0.25.multiExon.gtf" % (out_prefix),"w" )
      print "opening %s.FPKM0.5_rep0.25.multiExon.gtf" % (out_prefix)
      f_log = open( "select.log","w" )
      for gene in self.RPKMInfo['list_gene']:
         if gene not in self.RPKMInfo['gene_info']:
            print >>f_log, gene,1
            continue
         
         if 'is_pass' not in self.RPKMInfo['gene_info'][gene]:
            print >>f_log, gene,2
            continue
            
         if self.RPKMInfo['gene_info'][gene]['is_pass'] != 1:
            print >>f_log, gene,3
            continue
         
         if self.gene[ gene ]['is_intergenic'] == 0:
            print >>f_log, gene,4
            continue
         
         if len( self.gene[gene]['tran'] ) < 1:
            print >>f_log, gene,5
            continue
         
         for tids in self.gene[gene]['tran']:
            if len( self.gene[gene]['exon'][tids] ) > 1:
               for exon in self.gene[gene]['exon'][tids]:
                  if gene not in self.use_gene:
                     self.use_gene[gene] = 1
                  if tids not in self.use_tids:
                     self.use_tids[tids] = 1
                  print >>f_out, self.gene[gene]['line'][tids][exon]
      f_out.close()
      f_log.close()
      


   def get_gene_info(self):
      out_prefix = "%s"     % ( ".".join( self.gtf.split(".")[:-1] ) )
      print "Reading %s.genlen"                          % (out_prefix)
      f_in = open( "%s.genlen"                           % (out_prefix),"r" )
      f_out= open( "%s.FPKM0.5_rep0.25.multiExon.genlen" % (out_prefix),"w" )
      for line in f_in:
         line = line.strip('\n')
         f    = line.split()
         if f[1] in self.use_tids:
            print >>f_out, line
      f_in.close()
      f_out.close()
         
         
      print "Reading %s.gene.bed"                              % (out_prefix)
      f_in = open( "%s.gene.bed"                               % (out_prefix),"r" )
      f_out= open( "%s.FPKM0.5_rep0.25.multiExon.gtf.gene.bed" % (out_prefix),"w" )
      for line in f_in:
         line = line.strip('\n')
         f    = line.split()
         if f[3] in self.use_tids:
            print >>f_out, line
      f_in.close()
      f_out.close()
   
   
   def __load_tag(self):
      self.RPKMInfo['samp']['tag']  = []
      self.RPKMInfo['samp']['sam_type'] = {}
      self.RPKMInfo['samp']['sam_stage'] = {}
      self.RPKMInfo['samp']['sam_tissue'] = {}
      self.RPKMInfo['samp']['type_sam'] = {}
      self.RPKMInfo['samp']['type_stage_sam'] = {}
      self.RPKMInfo['samp']['stage_sam'] = {}
      self.RPKMInfo['samp']['tissue_sam'] = {}
      self.RPKMInfo['samp']['tissue_type_sam'] = {}
      self.RPKMInfo['samp']['tissue_type_stage_sam'] = {}
      self.RPKMInfo['samp']['tissue_stage_type_sam'] = {}
      for sam in self.RPKMInfo['samp']['list']:
         f      = sam.split("_")
         tissue = f[0]
         stage  = f[0]
         ltype  = "RNA"
         
         tag = "%s_%s" % (tissue,stage)
         if tag not in self.RPKMInfo['samp']['tag']:
            self.RPKMInfo['samp']['tag'].append(tag)

         self.RPKMInfo['samp']['sam_type'][sam]   = ltype
         self.RPKMInfo['samp']['sam_stage'][sam]  = stage
         self.RPKMInfo['samp']['sam_tissue'][sam] = tissue
         if ltype  not in self.RPKMInfo['samp']['type_sam']:
            self.RPKMInfo['samp']['type_sam'][ltype] = []
   
         if ltype  not in self.RPKMInfo['samp']['type_stage_sam']:
            self.RPKMInfo['samp']['type_stage_sam'][ltype] = {}
         if stage  not in self.RPKMInfo['samp']['type_stage_sam'][ltype]:
            self.RPKMInfo['samp']['type_stage_sam'][ltype][stage] = []
         
         
         if stage  not in self.RPKMInfo['samp']['stage_sam']:
            self.RPKMInfo['samp']['stage_sam'][stage] = []
         if tissue not in self.RPKMInfo['samp']['tissue_sam']:
            self.RPKMInfo['samp']['tissue_sam'][tissue] = []
         
         if tissue not in self.RPKMInfo['samp']['tissue_type_sam']:
            self.RPKMInfo['samp']['tissue_type_sam'][tissue] = {}
         if ltype  not in self.RPKMInfo['samp']['tissue_type_sam'][tissue]:
            self.RPKMInfo['samp']['tissue_type_sam'][tissue][ltype] = []
         
         if tissue not in self.RPKMInfo['samp']['tissue_type_stage_sam']:
            self.RPKMInfo['samp']['tissue_type_stage_sam'][tissue] = {}
         if ltype  not in self.RPKMInfo['samp']['tissue_type_stage_sam'][tissue]:
            self.RPKMInfo['samp']['tissue_type_stage_sam'][tissue][ltype] = {}
         if stage  not in self.RPKMInfo['samp']['tissue_type_stage_sam'][tissue][ltype]:
            self.RPKMInfo['samp']['tissue_type_stage_sam'][tissue][ltype][stage] = []
         
         if tissue not in self.RPKMInfo['samp']['tissue_stage_type_sam']:
            self.RPKMInfo['samp']['tissue_stage_type_sam'][tissue] = {}
         if stage  not in self.RPKMInfo['samp']['tissue_stage_type_sam'][tissue]:
            self.RPKMInfo['samp']['tissue_stage_type_sam'][tissue][stage] = {}
         if ltype  not in self.RPKMInfo['samp']['tissue_stage_type_sam'][tissue][stage]:
            self.RPKMInfo['samp']['tissue_stage_type_sam'][tissue][stage][ltype] = []
            
         self.RPKMInfo['samp']['type_sam'][ltype].append(sam)
         self.RPKMInfo['samp']['type_stage_sam'][ltype][stage].append(sam)
         self.RPKMInfo['samp']['stage_sam'][stage].append(sam)
         self.RPKMInfo['samp']['tissue_sam'][tissue].append(sam)
         self.RPKMInfo['samp']['tissue_type_sam'][tissue][ltype].append(sam)
         self.RPKMInfo['samp']['tissue_type_stage_sam'][tissue][ltype][stage].append(sam)
         self.RPKMInfo['samp']['tissue_stage_type_sam'][tissue][stage][ltype].append(sam)
      

   def __is_pass(self,line):
      """docstring for __is_pass"""
      M_tis_cnt = {}
      f = line.split()
      gene = f[0]
      
      is_pass = 0
      for tis in self.RPKMInfo['samp']['tissue_sam']:
         l_value = []
         for sam in self.RPKMInfo['samp']['tissue_sam'][tis]:
            sam_cnt = 0
            idx = self.RPKMInfo['samp']['list'].index(sam)
            value = float( f[idx+1] )
            l_value.append( value )
            if value < 0.25:
               sam_cnt += 1
         np_value = np.array( l_value )
         if np_value.mean() > 0.5 and sam_cnt == 0:
            is_pass = 1
      return is_pass



