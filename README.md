A pipeline which could processing from raw fastq reads to FPKM value and unique counts for each gene/repeat elements, RNA quantification using ERCC molecules and for basic statistics for mapping.

First, before this pipeline in a server, make sure the required modules were installed. If not, running the following scripts for deploying.
```bash
mkdir install_packages

### install python anaconda 2.2.0
cd software/
wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda-2.2.0-Linux-x86_64.sh
bash Anaconda-2.2.0-Linux-x86_64.sh  # prefix=/path/for/anaconda
mv Anaconda-2.2.0-Linux-x86_64.sh install_packages

### install R 3.2.0
wget http://cran.r-project.org/src/base/R-3/R-3.2.0.tar.gz
tar -zxvf R-3.2.0.tar.gz
cd R-3.2.0
./configure --prefix ~/software/R-3.2.0
make
make install
cd ..
mv R-3.2.0.tar.gz install_packages

### install samtools 0.1.18
### using old version because the latest one could have somewhat trouble with 
### other software like tophat.
wget http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2
tar -jxvf samtools-0.1.18.tar.bz2
cd samtools-0.1.18
make
cd ..
mv samtools-0.1.18.tar.bz2 install_packages

### install bwa 0.7.5a
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.5a.tar.bz2
tar -jxvf bwa-0.7.5a.tar.bz2
cd bwa-0.7.5a
make
cd ..
mv bwa-0.7.5a.tar.bz2 install_packages

### install bowtie2 2.2.3
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip
unzip bowtie2-2.2.3-linux-x86_64.zip
mv bowtie2-2.2.3-linux-x86_64.zip install_packages

### install tophat 2.0.12
wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.12.Linux_x86_64.tar.gz
tar -zxvf tophat-2.0.12.Linux_x86_64.tar.gz
mv tophat-2.0.12.Linux_x86_64.tar.gz install_packages

### install cufflinks 2.2.1
wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar -zxvf cufflinks-2.2.1.Linux_x86_64.tar.gz
mv cufflinks-2.2.1.Linux_x86_64.tar.gz install_packages

### install bedtools 2.24.0
git clone https://github.com/arq5x/bedtools2/
make

### install HTSeq
pip install HTSeq

### install tabix and bgzip
wget http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
tar -jxvf tabix-0.2.6.tar.bz2
cd tabix-0.2.6
make
cd ..
mv tabix-0.2.6.tar.bz2 install_packages

### install UCSC utilities
from http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/

```

After that, download this script:
```
cd $PYTHONPATH  # path for put the python packages. path/to/anaconda/lib/python2.7/site-packages/ for default
git clone https://github.com/hubqoaing/RNA_v2
```


Secondly, go to the ./setting file, and change the following values to your own path:
```python
self.Database       = "DIR/TO/DATABASE"          #line 56
self.sftw_py        = "DIR/TO/SOFTWARE_EXE_FILE" #line 78
self.sftw_pl        = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_tophat_dir= "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_cflk_dir  = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_bowtie_dir= "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_ucsc_dir  = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_samtools  = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_bedtools  = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_deseq     = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_bgzip     = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_tabix     = "DIR/TO/SOFTWARE_EXE_FILE"

```

Go to the analysis dictionary and copy the bin file here.
``` bash
cd PATH/FOR/ANALYSIS   # go to 
copy $PYTHONPATH/RNA_v2/run_mRNA.py ./
```

Next, make the input files. You can download these files in UCSC or so on and then using own-scripts to merge the ERCC information, and generate files in this format.
``` bash
vim sample_input.xls
==> sample_input.xls <==
sample		brief_name		stage			sample_group    ERCC_times	RFP_polyA	GFP_polyA	CRE_polyA	end_type	rename
NAME_FOR_RAW_FQ		NAME_FOR_PROCESSING	Group_FOR_STAGE	RNA             0.0    		0.0    		0.0		   0.0		 	PE          NAME_FOR_READING
```
Notice that only NAME_FOR_RAW_FQ were required that this NAME should be the same as 00.0.raw_fq/NAME.
NAME_FOR_PROCESSING will be the name for the rest analysis's results.
NAME_FOR_READING    will be the name for files in statinfo.
stage and sample_group could be writen as anything. It was here only for make the downstream analysis easily.

Before running this pipeline, put the fastq reads in the ./00.0.raw_data dictionary.
```bash
mkdir 00.0.raw_data
for i in `tail -n +2 sample_input.xls | awk '{print $1}`
do
    mkdir 00.0.raw_data/$i && ln -s PATH/TO/RAW_DATA/$i/*gz 00.0.raw_data/$i
done
```

After that, running this pipeline:
```bash
 python run_mRNA.py --ref YOUR_REF sample_input.xls
```

Wait for the results. 
Notice if you have to run it in a cluster, please do not running this scripts directly.
For example, if SGE system used, then:

Comments this command
```python
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
```
and using this command in modules in ./frame/*py
```
       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)
```
Method for submit jobs in other system were still developing.
