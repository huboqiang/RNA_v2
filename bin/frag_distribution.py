from __future__ import division
import re,sys,os,gzip

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pysam
import pandas as pd
import seaborn as sns
sns.set(style="ticks")

class refGeneLine(object):
    def __init__(self, f):
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


class FragDist(object):
    def __init__(self, infile_genelen, infile_refGene, infile_bam):
        
        self.infile_genlen = infile_genelen
        self.infile_refGene = infile_refGene
        self.infile_bam = infile_bam
        self.f_bam = pysam.Samfile(infile_bam,'rb')
        self.M_geneRead = {}

    def getLongestTran(self, infile_genlen):
        self.M_tran_used = {}
        with open(infile_genlen, "r") as f_genLen:
            for line in f_genLen:
                line = line.strip()
                f = line.split()
                if f[2] == f[3]:
                    self.M_tran_used[f[1]] = f[0]

    def read_refGene(self, infile_refGene, M_tran_used):
        with open(infile_refGene, "r") as f_refGene:
            for line in f_refGene:
                line = line.strip()
                f = line.split()
                if f[1] in self.M_tran_used:
                    self.M_geneRead[self.M_tran_used[f[1]]] =  self.getGeneReadInfo(f)
                    
    def getGeneReadInfo(self, f):
        m_refInfo = refGeneLine(f)
        record = self.f_bam.fetch(reference=m_refInfo.chrom, start=m_refInfo.txStart, end=m_refInfo.txEnd)
        for rec in record:
            read_begin = rec.pos
            read_end = rec.query_alignment_end