from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np



def get_gene_order(infile):
    f_infile = open(infile,"r")
    l_gene = []
    M_gene = {}
    for line in f_infile:
        line = line.strip('\n')
        f    = line.split()
        if f[0] not in M_gene:
            M_gene[f[0]] = 1
            l_gene.append(f[0])
    f_infile.close()
    return l_gene

def get_gene_FPKM(infile):
    f_infile = open(infile,"r")
    M_gene = {}
    f_infile.readline()
    for line in f_infile:
        line = line.strip('\n')
        f    = line.split()
        gene = f[0]
        FPKM = float(f[9])
        if gene not in M_gene:
            M_gene[gene] = []
        M_gene[gene].append(FPKM)
        
        gene = f[4]
        FPKM = float(f[9])
        if gene not in M_gene:
            M_gene[gene] = []
        M_gene[gene].append(FPKM)
        
    f_infile.close()
    return M_gene

def show_help():
    print >>sys.stderr,"\n\tpython",sys.argv[0],"/data/Analysis/huboqiang/tmp/prog_dev/RNA_v2/04.1.cufflinks/EPI_A3/genes.fpkm_tracking /data/Analysis/huboqiang/Database_RNA_v2/hg19/all.exon.sort.ERCC.gtf"
   


def main():
    try:
        infile = sys.argv[1]
        gen_gtf= sys.argv[2]
    except IndexError:
        show_help()
        sys.exit()
    
    
    gen_len = "%s.genlen" % ('.'.join(gen_gtf.split('.')[:-1]) )
    l_gene = get_gene_order(gen_len)
    M_gene = get_gene_FPKM(infile)
    
    for gene in l_gene:
        FPKM = 0
        if gene in M_gene:
            FPKM = sum(M_gene[gene])/len(M_gene[gene])
        
        print "%s\t%f" % (gene,FPKM)
        
if __name__ == '__main__':
    main()
