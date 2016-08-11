#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os

import cPickle as pickle
import numpy as np
from optparse   import OptionParser

import RNA_v2.frame.module01_mapping_from_raw   as m01
import RNA_v2.frame.module02_RNA_Quantification as m02
import RNA_v2.frame.module03_Stat_Sample        as m03

import RNA_v2.utils.module_create_database  as m_db

def prepare_optparser():
    usage ="""usage: %s [options] sample_input_information_xls

Reference fasta file and refGene file should be all put in dictionary
self.Database, which were defined in settings/projpath.py

If not put in self.Database, this program will downloading from UCSC.
Detail information could be reached in utils/module_create_database.py

Suport genome includes:
    http://hgdownload.soe.ucsc.edu/goldenPath


Using -h or --help for more information

Example:
    python %s  --ref hg19 samp_test.xls

    """ % (sys.argv[0],sys.argv[0])

    description = " mRNA analysis pipeline "

    optparser = OptionParser(
        version="%s v2.0 20150610" % (sys.argv[0]),
        description=description,
        usage=usage,
        add_help_option=False
    )

    optparser.add_option(
        "-r", "--ref", default="hg19",
        help="\nReference genome. [default: %default]"
    )

    optparser.add_option(
        "--extra_GTF", default=None,
        help="\nExtra GTF file for lncRNA. [default: %default]"
    )
    optparser.add_option(
        "--given_GTF", default=None,
        help="\nGiven GTF. Gene id for LncRNA should be lncGene"+\
        "while trans id should be lncTran. [default: %default]"
    )

    optparser.add_option(
        "-h", "--help", action="help",
        help="\nShow this help message and exit."
    )
    return optparser

def main():
    prepare_optparser()
    (options,args) = prepare_optparser().parse_args()
    try:
        sam_RNAinfo = args[0]
        ref       = options.ref
        extra_GTF = options.extra_GTF
        given_GTF = options.given_GTF
    except IndexError:
        prepare_optparser().print_help()
        sys.exit(1)

    part0 = m_db.DataBaseInit(ref, sam_RNAinfo, is_debug=0)
    part0.check_files()

    part1 = m01.Map_From_raw(ref, sam_RNAinfo, is_debug=0)
    part1.s01_QC(core_num=4)
    part1.s02_Tophat(core_num=1)

    part2 = m02.RNA_Quantification(ref, sam_RNAinfo, core_num=4, is_debug=0)
    part2.run_pipeline(extra_GTF, given_GTF, is_MergeSam=1)

    part3 = m03.SampStat(ref, sam_RNAinfo, given_GTF, is_debug=0)
    try:
        part3.Basic_Stat()
    except:
        pass
    try:
        part3.ERCC_Stat()
    except:
        pass
    part3.Repeat_Stat()

if __name__ == '__main__':
    main()
