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
from Crypto.Cipher import AES
import base64

import RNA_v2.utils.module_MongoDB as m_mgdb

f_in = open('/home/huboqiang/.ssh/id_rsa.pub.bak', "r")
id_rsa = f_in.readline().split()[1][0:32]
encoded = "TRoUpU6MPxd2csHTww1L6X+l5mpzzuy/JPMrvecUW1Pv0kGxXG3q11yxKhR0dVj2GnKgHvIwY2/XFOOMyGi1PwUB6ngWt/n9vXDsdUvgtrA="
f_in.close()

def check_gzipped(infile):
    f_infile = gzip.open(infile)
    gzipped = True
    try:
        f_infile.read( 10 )
    except IOError:
        gzipped = False
    f_infile.close()
    return gzipped

class MongoGTEx(object):
    """docstring for DbGene"""
    def __init__(self, infile_exp, infile_sam, head_skip=0, col_left=2, col_gene=2):
        self.infile_exp = infile_exp
        self.infile_sam = infile_sam
        is_gzipped = check_gzipped(self.infile_exp)
        if is_gzipped:
            self.f_exp = gzip.open(self.infile_exp, "r")
        else:
            self.f_exp = open(self.infile_exp, "r")
        
        for i in range(head_skip):
            self.f_exp.readline()
        
        self.pd_sam = pd.read_csv(infile_sam, sep="\t")
        self.mat_col_left = col_left
        self.mat_col_gene = col_gene
        
        cipher = AES.new(id_rsa, AES.MODE_ECB)
        uri = cipher.decrypt(base64.b64decode(encoded)).split()[0]
        self.client = MongoClient(uri)
        self.db_PGC_exp = self.client.myCollection_hg19_PGC_exp
        self.posts_PGC_exp = self.db_PGC_exp["GTEX_hg19"]
        self.__load_samInfo("SAMPID", "SMTS")

    def __load_samInfo(self, samID, tisID):
        self.struct = {}
        self.tis_sam = {}
        lineSam = self.f_exp.readline()
        lineSam = lineSam.strip()
        self.l_samExp = lineSam.split()[self.mat_col_left:]
        
        l_samples = list(self.pd_sam[samID].values)
        l_tissue = self.pd_sam[tisID].values
        M_samInfo = {}
        for i,sam in enumerate(self.l_samExp):
            idx_inExp = l_samples.index(sam)
            tissue = l_tissue[idx_inExp]
            if tissue is np.nan:
                continue
            if tissue not in self.struct:
                self.struct[tissue] = []
            if tissue not in self.tis_sam:
                self.tis_sam[tissue] = []
            self.struct[tissue].append(i)
            self.tis_sam[tissue].append(sam)
        
    def load_exp(self):
        for line in self.f_exp:
            line = line.strip()
            f = line.split()
            gene = f[self.mat_col_gene]
            INFO = {'gene': gene}
            for tissue in self.struct:
                np_idx = np.array(self.struct[tissue], dtype="int")
                np_val = np.array(f[self.mat_col_left:], dtype="float")
                tag = tissue
                INFO[tag] = list(np_val[np_idx])
            self.posts_PGC_exp.insert_one(INFO)
        self.posts_PGC_exp.create_index([("gene", 1)])
    
    def get_gene_pd(self, gene, l_tissue=None):
        post = self.db_PGC_exp['GTEX_hg19'].find_one({'gene':gene})
        if l_tissue is None:
            l_tissue = sorted(self.tis_sam.keys())
        
        l_sample = []
        l_stage = []
        l_val = []
        if post is not None:
            for tis in l_tissue:
                for i in range(len(self.tis_sam[tis])):
                    l_sample.append(self.tis_sam[tis][i])
                    l_stage.append(tis)
                    l_val.append(post[tis][i])

        pd_exp_plot = pd.DataFrame({
            "Sample": l_sample,
            "Stage" : l_stage,
            "exp"   : l_val
        })
        return pd_exp_plot
    
    def __del__(self):
        self.f_exp.close()
