from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np
from optparse   import OptionParser
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import scipy.stats
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch
import seaborn as sns
import pandas  as pd

import logging
logging.basicConfig(level=logging.INFO,format='%(levelname)-5s @ %(asctime)s: %(message)s ',stream=sys.stderr)


def get_stage_color( info_file ):
   pd_sampInfo = pd.read_csv( info_file,sep="\t",index_col=0 )
   l_stage     = []
   for stage in pd_sampInfo['stage']:
      if stage not in l_stage:
         l_stage.append( stage )
   
   stage_pal = range(0,len(l_stage))
   stage_lut = dict(zip(map(str, l_stage), stage_pal))
   stage_colors = np.array( pd_sampInfo['stage'].map(stage_lut) )
   return stage_colors
   
def get_tissue_color( info_file ):
   pd_sampInfo = pd.read_csv( info_file,sep="\t",index_col=0 )
   l_tissue     = []
   for tissue in pd_sampInfo['tissue']:
      if tissue not in l_tissue:
         l_tissue.append( tissue )
   
   tissue_pal = range(0,len(l_tissue))
   tissue_lut = dict(zip(map(str, l_tissue), tissue_pal))
   tissue_colors = np.array( pd_sampInfo['tissue'].map(tissue_lut) )
   
#   tissue_pal = sns.cubehelix_palette(len(l_tissue),light=.9, dark=.1, reverse=True, start=1, rot=-2)
#   tissue_lut = dict(zip(map(str, l_tissue), tissue_pal))
#   tissue_colors = np.array( pd_sampInfo['tissue'].map(tissue_lut) )
   return tissue_colors,l_tissue

   

def clean_axis(ax):
   """Remove ticks, tick labels, and frame from axis"""
   ax.get_xaxis().set_ticks([])
   ax.get_yaxis().set_ticks([])
   for sp in ax.spines.values():
      sp.set_visible(False)
        
class MatrixCorr(object):
   def __init__(self,infile,method="pearson"):
      self.infile = infile
      self.method = method
      
      self.out_cor_pdf   = "%s.corMat.%s.pdf" % ( ".".join( self.infile.split(".")[:-2] ), self.method )
      self.out_cor_xls   = "%s.corMat.%s.xls" % ( ".".join( self.infile.split(".")[:-2] ), self.method )
      
   def load_matrix(self):
      self.pd_mat = pd.read_csv( self.infile,sep="\t",index_col=[0] )
      print "A"
      self.pd_mat[ self.pd_mat<0.1 ] = 0.0
#      self.pd_mat = np.log10(self.pd_mat+1)
#      self.pd_mat.replace( 0.0,np.nan )
   
   def get_cor_matrix( self ):
      print "1"
      self.pd_cor   = self.pd_mat.corr( self.method )
      print "2"
      self.pd_cor.to_csv( self.out_cor_xls, sep="\t" )




   def plot_cor_matrix( self,info_file,cluster=1 ):
#      sns.set(style="darkgrid")
      self.pd_cor = pd.read_csv( self.out_cor_xls, sep="\t",index_col=0  )
      pd_sampInfo = pd.read_csv( info_file,sep="\t" )
      l_stage     = []
      for stage in pd_sampInfo['stage']:
         if stage not in l_stage:
            l_stage.append( stage )
      
      stage_pal = range(0,len(l_stage))
      stage_lut = dict(zip(map(str, l_stage), stage_pal))
      stage_colors = np.array( pd_sampInfo['stage'].map(stage_lut) )
      
      my_norm = matplotlib.colors.Normalize( 0,1 )
      
      
      dataMatrixOrdered = np.array(self.pd_cor.values)
      fig      = plt.figure(figsize=(15,15))
      
      
      
      if cluster == 1:
         # calculate pairwise distances for ROWS
         row_pairwise_dists = distance.squareform( distance.pdist(dataMatrixOrdered) )
         row_clusters = sch.linkage(row_pairwise_dists,method='ward')
         row_denAX = fig.add_axes([0.03,0.15,0.09,0.70])
         row_idx   = sch.dendrogram(row_clusters,color_threshold=np.inf,orientation='right')
         clean_axis(row_denAX)
         idx_row   = np.array(row_idx['leaves'])
         
         # calculate pairwise distances for COLUMNS
         col_pairwise_dists = distance.squareform( distance.pdist(dataMatrixOrdered.T) )
         col_clusters = sch.linkage(col_pairwise_dists,method='ward')
         col_denAX = fig.add_axes([0.20,0.86,0.70,0.09])
         col_idx   = sch.dendrogram(col_clusters,color_threshold=np.inf)
         clean_axis(col_denAX)
         idx_col   = np.array(col_idx['leaves'])
         
         dataMatrixOrdered_plot = dataMatrixOrdered[      idx_row,: ]
         dataMatrixOrdered_plot = dataMatrixOrdered_plot[ :,idx_col ]
         np_samp                = np.array(self.pd_cor.columns,dtype="string")
         np_samp_out            = np_samp[ idx_row ]
         stage_colors           = stage_colors[idx_col]
      else:
         dataMatrixOrdered_plot = dataMatrixOrdered
         np_samp_out            = np.array(self.pd_cor.columns,dtype="string")
         
         
      cmapMain = sns.diverging_palette(220, 10, s=90, l=30, as_cmap=True)
      
      cxmatrix = fig.add_axes([ 0.20, 0.15, 0.70, 0.70])
      cxmatrix.grid(True,color='black',linestyle='-',linewidth=3)
      
      cax02 = cxmatrix.imshow(dataMatrixOrdered_plot, aspect='auto', cmap=cmapMain,interpolation='nearest',norm=my_norm)
      cxmatrix.get_xaxis().set_ticks([])
      cxmatrix.get_yaxis().set_ticks([])
      cxmatrix.get_yaxis().set_ticks([])
      
      cxmatrix.get_xaxis().set_ticklabels([])
      for sp in cxmatrix.spines.values():
         sp.set_visible(False)

      for i in xrange( 0,len( np_samp_out ) ) :
         cxmatrix.text( len(np_samp_out),i,'%s' % ( np_samp_out[i]),size=7)

      
      cbaxes = fig.add_axes([0.03,0.90,0.10,0.02])
      cbar = plt.colorbar(cax02,cax=cbaxes, orientation='horizontal',ticks=[0,1])
      cbar.ax.yaxis.set_ticks_position('left')

      cbar.set_label("%s\nCorrelation" %  self.method, size=12)
      cbar.ax.tick_params(labelsize=10)
      
   
      ax1 = fig.add_axes([ 0.13, 0.15,  0.06, 0.7 ]  )
      ax1.get_xaxis().set_ticks( [] )
      ax1.get_yaxis().set_ticks( [] )
      ax1.get_yaxis().set_ticklabels( [] )
      ax1.get_xaxis().set_ticklabels( [] )
      ax1.tick_params(colors="white")
#      ax1.grid(True,color='white',linestyle='-',linewidth=3)

      np_show = np.array( stage_colors,dtype="float"  )
      cax1 = ax1.imshow( np.array([np_show]).T, aspect='auto', cmap='jet',interpolation='nearest')


      ax2 = fig.add_axes([ 0.20, 0.12,  0.70, 0.02 ]  )
      ax2.get_xaxis().set_ticks( [] )
      ax2.get_yaxis().set_ticks( [] )
      ax2.get_yaxis().set_ticklabels( [] )
      ax2.get_xaxis().set_ticklabels( [] )
      ax2.tick_params(colors="white")
      np_show = np.array( stage_colors,dtype="float"  )
      cax2 = ax2.imshow( np.array([np_show]), aspect='auto', cmap='jet',interpolation='nearest')
      
      for i in xrange( 0,len( np_samp_out ) ) :
         ax2.text( i+0.5,1,'%s' % ( np_samp_out[i]),rotation=270,ha = "right",size=7)
      
      
      pd_cor_new = pd.DataFrame( dataMatrixOrdered_plot, index=np_samp_out,columns=np_samp_out )

      file_prefix = ".".join( self.out_cor_xls.split(".")[:-1] )
      pd_cor_new.to_csv( "%s.cluster%d.xls" % (file_prefix,cluster), sep="\t" )

      plt.savefig("%s.cluster%d.pdf" % (file_prefix,cluster) ,format='pdf')
      
   
   def plot_PCA(self,info_file):
      from sklearn import decomposition
      from sklearn import datasets
      
      pd_sampInfo = pd.read_csv( info_file,sep="\t" )
      tissue_color,l_tissue = get_tissue_color( info_file )

      FPKM_mat = self.pd_mat.values.T
      pca = decomposition.PCA(n_components=6)
      pca.fit( FPKM_mat )
      np_pc = pca.transform(FPKM_mat)

      pd_PCA = pd.DataFrame( np_pc[:,0:2],index=self.pd_mat.columns,columns=["PC1","PC2"] )

      sns.set_style("white")
      fig = plt.figure(figsize=(30,15))
      ax = plt.subplot( 1, 1, 1 )
      
      
      tissue_pal = sns.color_palette("rainbow", len(l_tissue))

      tissue_lut = dict(zip(map(str, l_tissue), tissue_pal))
      
      plot_colors = np.array( pd_sampInfo['tissue'].map(tissue_lut) )
      for i,tis in enumerate(l_tissue):
         index = (tissue_color == i)
         ax.plot( np_pc[index,0],np_pc[index,1],".",color=tissue_pal[i],label=tis )
         
      for i,sam in enumerate( self.pd_mat.columns ):
         ax.text( np_pc[i,0],np_pc[i,1],"%s" % sam )
      
      sns.despine(trim=True)   
      
      box = ax.get_position()
      ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
      
      ax.legend(loc="center left",fontsize=9,bbox_to_anchor=(1, 0.5),fancybox=True, shadow=True)
      
      
      
      file_prefix = ".".join( self.out_cor_xls.split(".")[:-1] )
      
      pd_PCA.to_csv( "%s.sklearn.xls" % (file_prefix),sep="\t" )
      fig.savefig(   "%s.sklearn.pdf" % (file_prefix),format="pdf")
      
      
      
      
      
def prepare_optparser():
   usage ="""usage: %s [options]

   Using -h or --help for more information

Example:
   python %s merge.FPKM.test.less.xls  sample_name.20150509Reorder.less.xls  pearson 1

   """ % (sys.argv[0],sys.argv[0])

   description = "Printing ChIP seq correlation using read density (RPKM)."

   optparser = OptionParser(version="%s v0.1 20150131" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
   optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
   return optparser
   
def main():
   prepare_optparser()
   (options,args) = prepare_optparser().parse_args()
   try:
      infile    = args[0]
      info_file = args[1]
      method    = args[2]
      is_cluster= int(args[3])
   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)
   
   infile_dat = "%s.%s.dat" % ( ".".join( infile.split(".")[:-1] ), method )
   try:
      logging.info( 'Find %s.' % ( infile_dat ) )
      m_cor = pickle.load( open(infile_dat) )
      logging.info( 'Load %s done.' % ( infile_dat ) )
   except:
      logging.info( 'Did not find %s, trying to load matrix.' % ( infile_dat ) )
      m_cor = MatrixCorr( infile,method )
      m_cor.load_matrix()
      logging.info( 'Load matrix %s done. \nNow Calculating %s correlation.' % ( infile,method ) )
      m_cor.get_cor_matrix()
      logging.info( 'Calculating %s correlation for %s done.' % ( method,infile ) )
      pickle.dump( m_cor,open(infile_dat,"wb"),True )   
      logging.info( 'Write %s done.' % ( infile_dat ) )
   
   logging.info( 'Plotting %s %s correlation.' % ( infile,method ) )
   m_cor.plot_cor_matrix( info_file,is_cluster )

   logging.info( 'Plotting %s PCA.' % ( infile ) )
   m_cor.plot_PCA( info_file )

if __name__ == '__main__':
   main()