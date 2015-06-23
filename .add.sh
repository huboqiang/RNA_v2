git add  __init__.py
git commit -m "init"

git add run_mRNA.py
git commit -m "Command for running the whole pipeline."
 
git add .gitignore
git commit -m "ignore file types"

git add README.md
git commit -m "how to use the modules."
 
 
### submit bin dictionary ###
git add bin/multi-process.pl
git commit -m "multi-process shell running in nohup."
 
git add bin/qsub-sge.pl
git commit -m "multi-process shell running in sge-qsub system."

git add bin/s03_genePred2bed.py
git commit -m "Convert refGene(genePred) file to bed file."

git add bin/QC.pl
git commit -m "Quality control for reads before mapping."

git add bin/Repeat_Intersect2Count.py
git commit -m "Script used for repeat counting."

git add bin/grep_FPKM.py
git commit -m "Make FPKM value in cufflinks gene.track file in given order."

git add bin/find_ExonIntronIntergenic
git commit -m "Detect intron region, exon region and intergenic region. First define intron, then exon, rest for intergenic."

### submit database dictionary ###
git add database/ERCC.fa
git commit -m "Sequence of the ercc-spikein(92ercc+RGC)."
 
git add database/ERCC.gtf
git commit -m "Annotation for the ercc-spikein(92ercc+RGC)."

git add database/ercc.info.xls
git commit -m "Concerntration for the ercc-spikein(92ercc+RGC)."


### submit frame dictionary ###
git add frame/__init__.py
git commit -m "init."
 
git add frame/module01_mapping_from_raw.py
git commit -m "Infrastructure for mapping."

git add frame/module02_RNA_Quantification.py
git commit -m "Infrastructure for RNA quantification."

git add frame/module03_Stat_Sample.py
git commit -m "Infrastructure for Statistics of the final results."

### submit setting dictionary ###
git add settings/__init__.py
git commit -m "init."

git add settings/projpath.py
git commit -m "Setting for the path of the files for this project."

git add settings/scripts.py
git commit -m "Setting for how scripts were printed and executed using this pipeline."


### submit utilities module for this project ###
git add utils/__init__.py
git commit -m "init."

git add utils/module_CountInfo.py
git commit -m "module for processing with count matrix."

git add utils/module_create_database.py
git commit -m "Module for create database for this project."

git add utils/module_GTFFeature.py
git commit -m "Module for processing GTF file. Import for novo gene detection"

git add utils/module_Matrix.py
git commit -m "Module for FPKM-excel processing."

git add utils/module_refGene.py
git commit -m "Module for processing refGene genePred file."

git add utils/module_RepeatCount.py
git commit -m "Module for RPKM of repeat elements."

git add utils/module_running_jobs.py
git commit -m "module for jobs running."

git add utils/module_StatInfo.py
git commit -m "Model for dealing with statistics information on pipeline."


# 
# git add module01_mapping_from_sra.py
# git commit -m "SRA => fq => bam"
# 
# git add module_running_jobs.py

# 
# 
# git add module02_RNA_Quantification.py
# git commit -m "bam => excels for counts/FPKM."
# git add README.md
# git commit -m "Add guidance for running module02."

# git add module02_RNA_Quantification.py
# git commit -m "Fixed debug on reading Hisat log."
# git add module_RNAQuantPipe.py
# git commit -m "Fixed debug on reading Hisat log."
# git add README.md
# git commit -m "Add more files on file-system."
# 
# git add module_CountInfo.py
# git commit -m "Dealing with Count-matrix files preparing for DEG analysis, and calculate RPKM for novo-lincRNA to filter novo-lincRNAs with insufficient abundance."
# git add module_GTFFeature.py
# git commit -m "Reading and Parsing the GTF files."
# git add module_StatInfo.py
# git commit -m "Model for dealing with statistics information on pipeline."

#git add README.md
#git commit -m "Add more information."
#git add module01_mapping_from_sra.py
#git commit -m "Revise information on parameters for mapping."
#git add module02_RNA_Quantification.py
#git commit -m "Revise information on parameters for mapping."
#git add module_RNAQuantPipe.py
#git commit -m "Revise information on parameters for downstream analysis."
#git add module_StatInfo.py
#git commit -m "Revise module for sample information."


#git add bin/RNA01_plot_cor_density.py
#git commit -m "Plot RNA-FPKM corr-plot."
#git add bin/Repeat_Intersect2Count.py
#git commit -m "Count reads on repeats."
#git add bin/merge_2Dxls.py
#git commit -m "Merge FPKM-excels."
#git add bin/reorder_rev.py
#git commit -m "Reorder and remove samples in FPKM-excel."
#git add module03_Stat_Sample.py
#git commit -m "Module for mapping information."
#git add module_Matrix.py
#git commit -m "Module for FPKM-excel processing."
#git add module_RepeatCount.py
#git commit -m "Module for RPKM of repeat elements."

git push -u origin master
