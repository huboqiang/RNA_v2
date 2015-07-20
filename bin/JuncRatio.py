from __future__ import division
import re,sys,os
import cPickle as pickle
import subprocess,time
import numpy as np
from optparse   import OptionParser
import sys
import pysam

def is_intersect(beg, end, pos):
    idx = 0
    if pos == beg+1:
        idx = 1
    return idx
        
def prepare_optparser():
    usage ="""usage: %s [options] 

    Using -h or --help for more information

Example:
    python %s --ref_Splice /data/Analysis/huboqiang/Database_RNA_v2/mm9/refGene.all_splicing.sort.bed  --out_file /date/huboqiang/data_fgp_GSE44183_mouse/01.Tophat/mouse_2cell_1/accepted_hits.spliceRatio.bed /date/huboqiang/data_fgp_GSE44183_mouse/01.Tophat/mouse_2cell_1/accepted_hits.bam 
   
    """ % (sys.argv[0],sys.argv[0])

    description = " Get the Junction-reads Linear-reads ratio in junctions."

    optparser = OptionParser(version="%s v0.1 20150712" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
    optparser.add_option("-r", "--ref_Splice",    default="/data/Analysis/huboqiang/Database_RNA_v2/mm9/refGene.all_splicing.sort.bed",       help="\nPhreads uses, [default: %default]")
    optparser.add_option("-o", "--out_file", default="/data/Analysis/huboqiang/Database_RNA_v2/mm9/refGene.all_splicing.sort.bed",help="\nOutput file prefix[default: %default]")
    optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
    return optparser
    
def main():
    prepare_optparser()
    (options,args) = prepare_optparser().parse_args()
    try:
        raw_bam    = args[0]
        ref_Splice = options.ref_Splice
        out_file  = options.out_file
    except IndexError:
        prepare_optparser().print_help()
        sys.exit(1)
    
    if not os.path.isfile( "%s.bai" % (raw_bam) ):
        shell_info = "samtools index %s" % (raw_bam)
        print >>sys.stderr, shell_info
        p = subprocess.Popen(shell_info,shell='True')
        while 1:
            run_cnt = 0
            if p.poll() is None:
                run_cnt += 1
                time.sleep(3)
            if run_cnt == 0:
                break
    
    f_sam = pysam.Samfile( raw_bam,"rb" )
    f_refSplice = open( ref_Splice,"r" )
    f_out_file = open( out_file    ,"w" )
    
    total_circ   = 0
    pass_PE_circ = 0
    notP_PE_circ = 0
    
    total_read   = 0
    pass_PE_read = 0
    notP_PE_read = 0
    
    for line in f_refSplice:
        line  =  line.strip('\n')
        f   =  line.split()
        chrom =      f[0]
        beg   =  int(f[1])
        end   =  int(f[1])+1
                
        cnt_junc = 0
        cnt_linear = 0
        
        record = f_sam.fetch(reference=chrom, start=beg-1, end=beg+1)
        for rec in record:
            rec_beg = rec.reference_start
            idx     = rec_beg - beg
            for pair in rec.cigar:
                ctag = pair[0]
                leng = pair[1]
                pos  = beg+ idx
                
                is_linear = 0
                is_junction = 0
                if ctag == 0:
                    for i in xrange(leng):
                        is_overlap = is_intersect(beg, end, beg+idx)
                        if is_overlap:
                            is_linear = 1
                        idx += 1

                elif ctag == 1:
                    continue

                elif ctag == 2:
                    for i in xrange(leng):
                        idx += 1
                
                elif ctag == 3:
                    for i in xrange(leng):
                        is_overlap = is_intersect(beg, end, beg+idx)
                        if is_overlap:
#                            print rec, i
                            if i <2 or (leng-i)<2:
                                is_junction = 1
                        idx += 1
                
                cnt_junc += is_junction
                cnt_linear += is_linear
                
        
        print >>f_out_file, "%s\t%d\t%d" % (line, cnt_junc, cnt_linear)
            
            
    f_refSplice.close()
    f_sam.close()
              
    f_out_file.close()

if __name__ == '__main__':
    main()
