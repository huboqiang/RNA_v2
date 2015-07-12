from __future__ import division
import re,sys,os
import numpy as np
import subprocess
import time
from pybedtools import BedTool

def check_pos( chrom,beg,end,M_length ):
    if beg < 0:
        beg = 1
    if end > M_length[chrom]:
        end = M_length[chrom]-1
    return beg,end

def CpG_type(chrom,beg,end, in_fa,M_length):
   
    beg,end = check_pos(chrom,beg,end,M_length)
    seq = BedTool.seq((chrom, beg, end),in_fa)
    seq = seq.upper()
    i       = 0;
    status  = "ICP";
    maxRcpg = 0;
    while  i+500 < len(seq) :
        tmp = seq[i:i+500]
        C  = tmp.count("C")
        G  = tmp.count("G")
        CpG= tmp.count("CG")
        
        Cgc  = (C+G)/500
        Rcpg = 0
        if C != 0 and G != 0:
            Rcpg = 500*CpG/(C*G)
        
        if Rcpg >= 0.75 and Cgc >= 0.55:
            status = "HCP"
            break
        else:
            if maxRcpg < Rcpg:
                maxRcpg = Rcpg
        i += 5
   
    if maxRcpg < 0.48:
        status = "LCP"
    return status

def sort_bed(infile):
    shell_info = """
bedtools sort -i %s >%s.tmp && mv %s.tmp %s
    """ % (infile,infile,infile,infile)
    p = subprocess.Popen(shell_info,shell='True')
    while 1:
        run_cnt = 0
        if p.poll() is None:
            run_cnt += 1
            time.sleep(1)
        if run_cnt == 0:
            break

def gzip_file(infile):
    outfile = ".".join(infile.split(".")[:-1])
    shell_info = """
zcat %s >%s
    """ % (infile, outfile)
    p = subprocess.Popen(shell_info,shell='True')
    while 1:
        run_cnt = 0
        if p.poll() is None:
            run_cnt += 1
            time.sleep(1)
        if run_cnt == 0:
            break
    
    return outfile

def put_pos(info, line):
    f = line.split()
    info['chr'] = f[0]
    info['beg'] = f[1]
    info['end'] = f[2]
    info['str'].append(f[3])
    info['gen'].append(f[4])
    info['tran'].append(f[5])
    return info

def report_inf(info):
    l_out = [
        info['chr'],
        info['beg'],
        info['end'],
        ",".join(list(set(info['str']))),
        ",".join(list(set(info['gen']))),
        ",".join(list(set(info['tran']))),
    ]
    return "\t".join(l_out)

def clean_info():
    info = {
        'chr' : "",
        'beg' : -1,
        'end' : -1,
        'str' : [],
        'gen' : [],
        'tran': []
    }
    return info

def uniq_bed(infile, col=5):
    in_prefix = ".".join( infile.split(".")[:-1] )
    out_file  = "%s.uniq_pos.bed" % in_prefix
    f_infile  = open(infile  , "r")
    f_outfile = open(out_file, "w")
    info = clean_info()
    pre_inf = ""
    for line in f_infile:
        line = line.strip('\n')
        f    = line.split()
        
        pos_inf = "_".join(f[0:3])
        if pos_inf == pre_inf:
            info = put_pos(info, line)
        else:
            if pre_inf != "":
                out = report_inf(info)
                print >>f_outfile, out
                info = clean_info()
            info = put_pos(info, line)
        pre_inf = pos_inf

    out = report_inf(info)
    print >>f_outfile, out
    f_outfile.close()            
                

class refGeneLine(object):
    def __init__(self,line):
        f = line.split()
        self.bin	      = int(f[0])
        self.name	      =     f[1]
        self.chrom	      =     f[2]
        self.strand       =     f[3]
        self.txStart      = int(f[4])
        self.txEnd        = int(f[5])
        self.cdsStart     = int(f[6])
        self.cdsEnd       = int(f[7])
        self.exonCount    = int(f[8])
        self.exonStarts   = np.array( f[ 9].split(",")[:-1],dtype=int )
        self.exonEnds     = np.array( f[10].split(",")[:-1],dtype=int )
        self.score        = int(f[11])
        self.name2        =     f[12]
        self.cdsStartStat =     f[13]
        self.cdsEndStat   =     f[14]
        self.exonFrames   =     f[15]

class RefGeneTxt(object):
    def __init__(self,infile,genomeFa):
        self.infile   = infile
        self.genomeFa = genomeFa
        self.length   = {}
        self.__load_length()
        
        if self.infile.split(".")[-1] == "gz":
            gzipped_file = gzip_file(self.infile)
            self.infile = gzipped_file
        
        self.prefix = ".".join( self.infile.split(".")[:-1] )
        self.infile_sort = "%s.sort.txt" % ( self.prefix )
        if not os.path.isfile( self.infile_sort ):
            self.__sort_refgene()
      
    def __sort_refgene(self):
        shell_info = """
for i in chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19\\
    chr2 chr20 chr21 chr22 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrX chrY 
do
    sort -k5n %s | grep -w $i
done >%s
        """ % ( self.infile, self.infile_sort )
        p = subprocess.Popen(shell_info,shell='True')
        while 1:
            run_cnt = 0
            if p.poll() is None:
                run_cnt += 1
                time.sleep(1)
            if run_cnt == 0:
                break
            
    def __load_length(self):
        f_inFai = open( "%s.fai" % (self.genomeFa),"r")
        for line in f_inFai:
            line = line.strip('\n')
            f    = line.split()
            self.length[ f[0] ] = int( f[1] )
        f_inFai.close()
    
    def RefSplicing(self):
        all_splice_site    = "%s.all_splicing.xls"    % (self.prefix)        
        f_infile_sort     =  open(self.infile_sort,"r")
        
        f_all_splice_site    = open(all_splice_site    ,"w")
        
        self.gen_splice = {}
        for line in f_infile_sort:
            line = line.strip('\n')
            line_info = refGeneLine(line)
            chrom = line_info.chrom
            gen   = line_info.name2
            tid   = line_info.name
            np_beg = line_info.exonStarts
            np_end = line_info.exonEnds
            
            if gen not in self.gen_splice:
                self.gen_splice[gen] = {}
            
            self.gen_splice[gen][tid] = {
                'chr': chrom, 
                'beg': np_beg, 
                'end': np_end,
                'str': line_info.strand,
            }
        f_infile_sort.close()
        
        for gen in sorted(self.gen_splice):
            tid_len = len(self.gen_splice[gen])
            l_pos_cnt  = {}
            l_pos_type = {}
            l_tid = self.gen_splice[gen].keys()
            for tid in self.gen_splice[gen]:
                for pos in self.gen_splice[gen][tid]['beg']:
                    if pos not in l_pos_cnt:
                        l_pos_cnt[pos]  = 0
                        l_pos_type[pos] = 'beg'
                        if line_info.strand == "-":
                            l_pos_type[pos] = 'end'
                    l_pos_cnt[pos] += 1
                    
                for pos in self.gen_splice[gen][tid]['end']:
                    if pos not in l_pos_cnt:
                        l_pos_cnt[pos] = 0
                        l_pos_type[pos] = 'end'
                        if line_info.strand == "-":
                            l_pos_type[pos] = 'beg'
                    l_pos_cnt[pos] += 1
                    
            for pos in l_pos_cnt:
                chrom  = self.gen_splice[gen][l_tid[0]]['chr']
                strand = self.gen_splice[gen][l_tid[0]]['str']
                
                stat = "diff"
                if l_pos_cnt[pos] == tid_len:
                    ### this splicing sites appeared in all samples. Common
                    stat = "common"
                
                l_out = [
                    chrom,
                    pos,
                    pos,
                    strand,
                    l_pos_type[pos],
                    stat,
                    gen
                ]
                l_out = [ str(s) for s in l_out ]
                print >>f_all_splice_site, "\t".join(l_out)
        
        f_all_splice_site.close()
    
    def refGene2bed(self,up_stream,down_stream,ext_type=""):
        self.refGene_ext_bed = "%s.up%d_down%d.%sBsorted.bed"               %\
            (self.prefix, up_stream, down_stream, ext_type)
            
        f_infile_sort     =  open(self.infile_sort,"r")
        f_refGene_ext_bed =  open(self.refGene_ext_bed,"w")
        
        for line in f_infile_sort:
            line = line.strip('\n')
            line_info = refGeneLine(line)
            chrom = line_info.chrom
            
            beg = line_info.txStart - up_stream
            end = line_info.txEnd   + down_stream
            if line_info.strand == "-":
                beg = line_info.txStart   - down_stream
                end = line_info.txEnd     + up_stream
            
            if ext_type == "tss." or ext_type == "promoter.":
                beg = line_info.txStart - up_stream
                end = line_info.txStart + down_stream
                if line_info.strand == "-":
                    beg = line_info.txEnd - down_stream
                    end = line_info.txEnd + up_stream
              
           
            strand = line_info.strand
            tid    = line_info.name
            ltype  = "protein_coding"
            if tid[0:2] == "NR":
                ltype  = "noncoding"
            
            gid = line_info.name2
            out = "%s\t%d\t%d\t%s\t%s\t%s\t%s"                              %\
                (chrom,beg,end, strand,tid,ltype,gid)
               
            print >>f_refGene_ext_bed, out
        f_infile_sort.close()
        f_refGene_ext_bed.close()
        
    def get_junc(self, ext_l=1):
        self.refGene_ext_bed = "%s.junc.ext_len%d.bed" % (self.prefix, ext_l)
            
        f_infile_sort     =  open(self.infile_sort,"r")
        f_refGene_ext_bed =  open(self.refGene_ext_bed,"w")
        
        for line in f_infile_sort:
            line = line.strip('\n')
            line_info = refGeneLine(line)
            chrom = line_info.chrom
            
            l_beg = line_info.exonStarts
            l_end = line_info.exonEnds
            
            if len(l_beg)>1:
                if line_info.strand == "+":
                    for i in xrange(len(l_beg)-1):
                        exon  = i + 1
                        l_out1 = [
                            chrom,
                            l_end[i] - ext_l,
                            l_end[i] + ext_l,
                            line_info.strand,
                            line_info.name2,
                            "%s_junc%d_l" % (line_info.name, exon)
                        ]
                        l_out2 = [
                            chrom,
                            l_beg[i+1] - ext_l,
                            l_beg[i+1] + ext_l,
                            line_info.strand,
                            line_info.name2,
                            "%s_junc%d_r" % (line_info.name, exon)
                        ]
                        l_out1 = [ str(i) for i in l_out1 ]
                        l_out2 = [ str(i) for i in l_out2 ]
                        print >>f_refGene_ext_bed, "\t".join(l_out1)
                        print >>f_refGene_ext_bed, "\t".join(l_out2)
                else:
                    for i in xrange(len(l_beg)-1):
                        exon = len(l_beg) - 1 - i
                        l_out1 = [
                            chrom,
                            l_end[i] - ext_l,
                            l_end[i] + ext_l,
                            line_info.strand,
                            line_info.name2,
                            "%s_junc%d_r" % (line_info.name, exon)
                        ]
                        l_out2 = [
                            chrom,
                            l_beg[i+1] - ext_l,
                            l_beg[i+1] + ext_l,
                            line_info.strand,
                            line_info.name2,
                            "%s_junc%d_l" % (line_info.name, exon),
                        ]
                        l_out1 = [ str(i) for i in l_out1 ]
                        l_out2 = [ str(i) for i in l_out2 ]
                        print >>f_refGene_ext_bed, "\t".join(l_out1)
                        print >>f_refGene_ext_bed, "\t".join(l_out2)
                
        f_infile_sort.close()
        f_refGene_ext_bed.close()

        sort_bed(self.refGene_ext_bed)
        uniq_bed(self.refGene_ext_bed, col=5)
        
   
    def refGeneInfo(self,TSS_up_stream=1000,TSS_down_stream=500):
        self.tranInfo = {}
        self.geneInfo = {}
      
        f_infile_sort = open( self.infile_sort,"r" )
      
        for line in f_infile_sort:
            line = line.strip('\n')
            line_info = refGeneLine( line )
            chrom = line_info.chrom
            beg   = line_info.txStart
            end   = line_info.txEnd
            
            strand = line_info.strand
            tid    = line_info.name
            ltype  = "protein_coding"
            if tid[0:2] == "NR":
                ltype  = "noncoding"

            gid    = line_info.name2
            
            tss_beg     = line_info.txStart - TSS_up_stream
            tss_end     = line_info.txStart + TSS_down_stream
            if line_info.strand == "-":
                tss_beg = line_info.txEnd   - TSS_down_stream
                tss_end = line_info.txEnd   + TSS_up_stream
               
#            tid_type = CpG_type(chrom, tss_beg, tss_end, 
#                                self.genomeFa, self.length)
            
            if gid not in self.geneInfo:
                self.geneInfo[ gid ] = { 
                    'tran':[],  'line':{},  'max_len':0, 
                    'max_len_tid':"None",'max_len_tid_type':""
                }
            
            self.geneInfo[ gid ][ 'tran' ].append( tid )
            self.geneInfo[ gid ][ 'line' ][ tid ] = line
            if end - beg > self.geneInfo[ gid ][ 'max_len' ]:
                self.geneInfo[ gid ][ 'max_len'     ] = end-beg
                self.geneInfo[ gid ][ 'max_len_tid' ] = tid
                self.geneInfo[ gid ][ 'max_len_tid_type' ] = tid_type
               
            if tid not in self.tranInfo:
                self.tranInfo[ tid ] = { 
                    'beg':beg, 'end':end,'chr':chrom,
                    'gene':gid, 'tmp_len':end - beg, 
                    'strand':strand,'tid_type':tid_type,
                    'exon_cnt':len(line_info.exonStarts), 
                    'exon_pos':{ 'left':line_info.exonStarts,
                                'right':line_info.exonEnds,
                               },
                }
        f_infile_sort.close()
            
    def Only_LongestTid_Bed(self,up_stream,down_stream,ext_type=""):
        infile = "%s.up%d_down%d.%sBsorted.bed"                             %\
                 (self.prefix, up_stream, down_stream, ext_type)
                 
        otfile = "%s.up%d_down%d.%sBsorted.longestTid.bed"                  %\
                 (self.prefix, up_stream, down_stream, ext_type)
                 
        f_infile = open(infile,"r")
        f_otfile = open(otfile,"w")
        for line in f_infile:
            line = line.strip('\n')
            f    = line.split()
            tid  = f[4]
            gid  = f[6]
            if tid == self.geneInfo[ gid ][ 'max_len_tid' ]:
                print >>f_otfile, "%s\t%s"                                  %\
                                    (line,self.tranInfo[tid]['tid_type'])

        f_infile.close()
        f_otfile.close()
        sort_bed( otfile )
      
