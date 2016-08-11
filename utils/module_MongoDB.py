#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import gzip
import subprocess
import time
import string

import numpy as np
import pandas as pd
from scipy import stats
import functools
from pymongo import MongoClient
from bson.objectid import ObjectId
from optparse import OptionParser
import logging
import tabix

logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
    stream=sys.stderr
)

class DbGene(object):
    """docstring for DbGene"""
    def __init__(self):
        super(DbGene, self).__init__()

    def loadExp_from_file(self, infile_exp, infile_sam):
        self.f_exp = gzip.open(infile_exp, "r")
        self.pd_sam = pd.read_csv(infile_sam, sep="\t")
        self.__load_samInfo()
        self.__load_exp()
    
    def __load_samInfo(self):
        self.struct = {}
        self.f_exp.readline()
        self.f_exp.readline()
        lineSam = self.f_exp.readline()
        lineSam = lineSam.strip()
        l_samExp = lineSam.split()[2:]
        
        l_samples = list(self.pd_sam['SAMPID'].values)
        l_tissue = self.pd_sam['SMTS'].values
        l_samExpIdx = []
        M_samInfo = {}
        for i,sam in enumerate(l_samExp):
#            print sam, l_samExp
            idx_inExp = l_samples.index(sam)
            tissue = l_tissue[idx_inExp]
            l_samExpIdx.append(i)
            if tissue is np.nan:
                continue
            if tissue not in self.struct:
                self.struct[tissue] = []

            self.struct[tissue].append(i)
    
    def __load_exp(self):
        for line in self.f_exp:
            line = line.strip()
            f = line.split()
            gene = f[1]
            INFO = {'gene': gene}
            for tissue in self.struct:
                np_idx = np.array(self.struct[tissue], dtype="int")
                np_val = np.array(f[2:], dtype="float")
                tag = tissue
                INFO[tag] = list(np_val[np_idx])
            self.posts_PGC_exp.insert_one(INFO)

        self.posts_PGC_exp.create_index([("gene", 1)])
