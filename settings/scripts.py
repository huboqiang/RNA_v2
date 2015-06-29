#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import subprocess
import cPickle as pickle
import time

import RNA_v2.utils.module_running_jobs as m_jobs
import RNA_v2.settings.projpath         as m_proj

class Scripts(m_proj.ProjInfo):
    def __init__(self):
        super( Scripts,self ).__init__()
        
    def define_files(self,ref):
        dir_db = "%s/%s" % (self.Database, ref)
        self.ref           = ref
        self.bt_index_base = "%s/%s_ERCC92_RGC"         % (dir_db, ref)
        self.genome_gtf    = "%s/refGene.gtf"           % (dir_db)
        self.intragenic_bed= "%s/region.Intragenic.bed" % (dir_db)
        self.rmsk_bed      = "%s/chrom.sort.bed"        % (dir_db)
        self.ercc_info_file= "%s/database/ercc.info.xls"% (self.path)

    def db_01_DownloadRef(self):
        l_sh_info = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("dir_path=%s"          % (self.path))
        l_sh_info.append("""
cd $dir_database

wget http://hgdownload.soe.ucsc.edu/goldenPath/${ref}/bigZips/chromFa.tar.gz

tar -zxvf $dir_database/chromFa.tar.gz

for i in {1..22} X Y M
do
    cat $dir_database/chr$i.fa
done  >$dir_database/${ref}_ERCC92_RGC.fa && rm $dir_database/chr*fa      &&\\

cat $dir_path/database/ERCC.fa >>$dir_database/${ref}_ERCC92_RGC.fa
        """)
        return l_sh_info

    def db_02_BuildRefIndex(self):
        l_sh_info = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("bowtie_dir=%s"   % (self.sftw_bowtie_dir))
        l_sh_info.append("samtools_exe=%s" % (self.sftw_samtools))
        l_sh_info.append("""
bowtie_build_exe=$bowtie_dir/bowtie2-build

$samtools_exe faidx $dir_database/${ref}_ERCC92_RGC.fa

$bowtie_build_exe $dir_database/${ref}_ERCC92_RGC.fa                        \\
    $dir_database/${ref}_ERCC92_RGC
        """)
        return l_sh_info
    
    def db_03_RefGene(self):
        l_sh_info = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("bedtools_exe=%s" % (self.sftw_bedtools))
        l_sh_info.append("ucsc_dir=%s"     % (self.sftw_ucsc_dir))
        l_sh_info.append("bin=%s"          % (self.bin))
        l_sh_info.append("dir_path=%s"     % (self.path))
        l_sh_info.append("""
cd $dir_database
wget http://hgdownload.soe.ucsc.edu/goldenPath/${ref}/database/refGene.txt.gz

### remove chromosome fragments(unassembled).
for i in {1..22} X Y M
do
    zcat $dir_database/refGene.txt.gz | grep -w chr$i
done >$dir_database/tmp
gzip $dir_database/tmp
mv $dir_database/tmp.gz $dir_database/refGene.txt.gz

# ref.sort.txt
# For CIRCExplorer(circular RNA pipeline)
zcat $dir_database/refGene.txt.gz                                          |\\
awk '{
    OFS="\\t";
    print $13,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
}' /dev/stdin                                                              |\\
sort /dev/stdin >$dir_database/ref.sort.txt                               &&\\

# refGene.bed
zcat $dir_database/refGene.txt.gz                                          |\\
awk '{
    tag="noncoding";
    if($4~/^NM/){tag="protein_coding"};
    OFS="\\t";
    print $3,$5,$6,$2,$4,$10,$11,tag,$13
}' /dev/stdin                                                              |\\
python $bin/s03_genePred2bed.py /dev/stdin                                 |\\
$bedtools_exe sort -i /dev/stdin >$dir_database/refGene.bed               &&\\

# region.Intragenic.bed
# For novo lncRNA detection
$bin/find_ExonIntronIntergenic/find_ExonIntronIntergenic                    \\
    $dir_database/refGene.bed                                               \\
    $dir_database/${ref}_ERCC92_RGC.fa.fai >$dir_database/pos.bed         &&\\

grep -v "Intergenic" $dir_database/pos.bed                                 |\\
    awk '{OFS="\t";print $1,$2,$3,"Intragenic"}' /dev/stdin                 \\
    >$dir_database/region.Intragenic.bed                                  &&\\

# refGene.gtf
# For mapping
zcat $dir_database/refGene.txt.gz                                          |\\
cut -f 2-                                                                  |\\
$ucsc_dir/genePredToGtf file stdin /dev/stdout                             |\\
grep -w exon                                                               |\\
$bedtools_exe sort -i /dev/stdin >$dir_database/refGene.gtf               &&\\
cat $dir_path/database/ERCC.gtf >>$dir_database/refGene.gtf
        """)
        return l_sh_info

    def db_04_rmsk(self):
        l_sh_info = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("bedtools_exe=%s" % (self.sftw_bedtools))
        l_sh_info.append("ucsc_dir=%s"     % (self.sftw_ucsc_dir))
        l_sh_info.append("bin=%s"          % (self.bin))
        l_sh_info.append("dir_path=%s"     % (self.path))
        l_sh_info.append("""
cd $dir_database
wget http://hgdownload.soe.ucsc.edu/goldenPath/${ref}/database/rmsk.txt.gz

zcat $dir_database/rmsk.txt.gz                                             |\\
awk '{
    OFS="\\t";
    print $6,$7,$8,$2,".",".",".","("$9")",$10,$11,$12 "/" $13,$14,$15,$16,$17
}' /dev/stdin                                                              |\\
tail -n +2  /dev/stdin >$dir_database/chrom.bed

for i in {1..22} X Y M
do
    grep -w chr$i $dir_database/chrom.bed
done >$dir_database/tmp
mv $dir_database/tmp $dir_database/chrom.bed

$bedtools_exe sort -i   $dir_database/chrom.bed >$dir_database/chrom.sort.bed
        """)
        return l_sh_info


    def db_05_Transcriptome(self):
        l_sh_info = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("tophat_dir=%s"   % (self.sftw_tophat_dir))
        l_sh_info.append("bowtie_dir=%s"   % (self.sftw_bowtie_dir))
        l_sh_info.append("""
bowtie_build_exe=$bowtie_dir/bowtie2-build
$tophat_dir/gtf_to_fasta $dir_database/refGene.gtf                          \\
    $dir_database/${ref}_ERCC92_RGC.fa                                      \\
    $dir_database/${ref}_ERCC92_RGC.refGene.fa

$bowtie_build_exe $dir_database/${ref}_ERCC92_RGC.refGene.fa                \\
    $dir_database/${ref}_ERCC92_RGC.refGene
        """)
        return l_sh_info




    def s_01_QC(self):
        l_sh_info     = []
        l_sh_info.append("samp=$1")
        l_sh_info.append("data_type=$2")
        l_sh_info.append("pl_exe=%s"       % (self.sftw_pl) )
        l_sh_info.append("pl_QC=%s/QC.pl"  % (self.bin)     )
        l_sh_info.append("in_dir=%s"       % (self.dir_raw_data  ) )
        l_sh_info.append("out_dir=%s"      % (self.dir_clean_data) )
        l_sh_info.append("""""")
        l_sh_info.append("""
$pl_exe $pl_QC --indir $in_dir --outdir $out_dir --sample $samp             \\
    --end $data_type""")
    
        return l_sh_info
        
    def s_02_Tophat(self):
        l_sh_info    = []
        l_sh_info.append("samp_name=$1")
        l_sh_info.append("brief_name=$2")
        l_sh_info.append("data_dype=$3")
        l_sh_info.append("tophat_py=%s/tophat"  % (self.sftw_tophat_dir))
        l_sh_info.append("cln_dir=%s"           % (self.dir_clean_data))
        l_sh_info.append("tophat_dir=%s"        % (self.dir_tophat))
        l_sh_info.append("genome=%s"            % (self.bt_index_base))
        l_sh_info.append("gtf_file=%s"          % (self.genome_gtf))
        l_sh_info.append("""
if [ $data_dype == "1" ]
    then $tophat_py                                                         \\
           -p 8 -G $gtf_file --library-type fr-unstranded                   \\
           --transcriptome-index $genome.refGene                            \\
           -o $tophat_dir/$brief_name  $genome                              \\
           $cln_dir/$samp_name/1.cln.fq.gz
fi
if [ $data_dype == "2" ]
    then $tophat_py                                                         \\
           -p 8 -G $gtf_file --library-type fr-unstranded                   \\
           --transcriptome-index $genome.refGene                            \\
           -o $tophat_dir/$brief_name  $genome                              \\
           $cln_dir/$samp_name/1.cln.fq.gz $cln_dir/$samp_name/2.cln.fq.gz
fi
        """)
        return l_sh_info


    def s_03_HTSeq_known(self):
        l_sh_info    = []
        l_sh_info.append("samp_name=$1")
        l_sh_info.append("py_exe=%s"        % (self.sftw_py))
        l_sh_info.append("samtools_exe=%s"  % (self.sftw_samtools))
        l_sh_info.append("deseq_exe=%s"     % (self.sftw_deseq))
        l_sh_info.append("tophat_dir=%s"    % (self.dir_tophat))
        l_sh_info.append("HTS_k_dir=%s"     % (self.dir_HTS_known))
        l_sh_info.append("known_GTF=%s"     % (self.genome_gtf))
        l_sh_info.append("""
$samtools_exe view  -H                                                      \\
    $tophat_dir/$samp_name/accepted_hits.bam                                \\
   >$tophat_dir/$samp_name/accepted_hits.header.sam

$samtools_exe sort  -n -m 200000000                                         \\
    $tophat_dir/$samp_name/accepted_hits.bam                                \\
    $tophat_dir/$samp_name/accepted_hits.sort_name

$samtools_exe view  -o                                                      \\
    $tophat_dir/$samp_name/accepted_hits.sort_name.sam                      \\
    $tophat_dir/$samp_name/accepted_hits.sort_name.bam 

[ ! -d $HTS_k_dir/$samp_name ] && mkdir -p $HTS_k_dir/$samp_name

$py_exe $deseq_exe                                                          \\
    -s no -f sam -a 10                                                      \\
    -o $tophat_dir/$samp_name/accepted_hits.sort_name.gene.sam              \\
    $tophat_dir/$samp_name/accepted_hits.sort_name.sam  $known_GTF          \\
    >$HTS_k_dir/$samp_name/$samp_name.dexseq.txt                          &&\\

grep -v -P '^ERCC-|^RGC-|MIR|SNORD|Mir|Snord'                               \\
    $HTS_k_dir/$samp_name/$samp_name.dexseq.txt                             \\
    >$HTS_k_dir/$samp_name/$samp_name.dexseq_clean.txt			          &&\\
    
grep -P '^ERCC-|^RGC-'                                                      \\
    $HTS_k_dir/$samp_name/$samp_name.dexseq.txt                             \\
    >$HTS_k_dir/$samp_name/$samp_name.dexseq_ERCC_RGCPloyA.txt	          &&\\

grep "__no_feature"                                                         \\
    $tophat_dir/$samp_name/accepted_hits.sort_name.gene.sam                |\\
    grep -v chrM                                                           |\\
	cat  $tophat_dir/$samp_name/accepted_hits.header.sam /dev/stdin        |\\
	$samtools_exe view -Sb /dev/stdin                                       \\
    >$tophat_dir/$samp_name/accepted_hits.genome.bam                      &&\\
    
$samtools_exe sort  -m 200000000                                            \\
    $tophat_dir/$samp_name/accepted_hits.genome.bam                         \\
    $tophat_dir/$samp_name/accepted_hits.genome.sort                      &&\\
    
rm  $tophat_dir/$samp_name/accepted_hits.sort_name.sam                      \\
    $tophat_dir/$samp_name/accepted_hits.sort_name.bam                      \\
    $tophat_dir/$samp_name/accepted_hits.sort_name.gene.sam                 \\
    $tophat_dir/$samp_name/accepted_hits.genome.bam
        """)
        return l_sh_info
        
        
        
    def s_04_1_novo_HTSeq_unknown(self):
        l_sh_info    = []
        l_sh_info.append("brief_name=$1")
        l_sh_info.append("sftw_cflk_dir=%s"    % (self.sftw_cflk_dir))
        l_sh_info.append("tophat_dir=%s"       % (self.dir_tophat))
        l_sh_info.append("dir_cufflinks_uk=%s" % (self.dir_cufflinks_unknown))
        l_sh_info.append("""
in_bam=$tophat_dir/$brief_name/accepted_hits.genome.sort.bam
out_dir=$dir_cufflinks_uk/$brief_name

$sftw_cflk_dir/cufflinks                                                    \\
   -p 8  -u                                                                 \\
   -o $out_dir                                                              \\
   $in_bam
        """)
        return l_sh_info
        
        
    def s_04_2_novo_cuffcomp(self):
        l_sh_info    = []
        l_sh_info.append("dir_cufflinks_uk=%s" % (self.dir_cufflinks_unknown))
        l_sh_info.append("sftw_cflk_dir=%s"    % (self.sftw_cflk_dir))
        l_sh_info.append("""
out_prefix=$1
shift

$sftw_cflk_dir/cuffcompare    \\
   -o  $out_prefix       \\
   -T  $@                \\
        """)
        return l_sh_info
        
        
    def s_04_3_novo_HTSeq_unknown(self):
        l_sh_info    = []
        l_sh_info.append("samp_name=$1")
        l_sh_info.append("unknown_GTF=$2")
        l_sh_info.append("py_exe=%s"           % (self.sftw_py))
        l_sh_info.append("samtools_exe=%s"     % (self.sftw_samtools))
        l_sh_info.append("deseq_exe=%s"        % (self.sftw_deseq))
        l_sh_info.append("tophat_dir=%s"       % (self.dir_tophat))
        l_sh_info.append("HTS_u_dir=%s"        % (self.dir_HTS_unknown))
        l_sh_info.append("sftw_cflk_dir=%s"    % (self.sftw_cflk_dir))
        l_sh_info.append("tophat_dir=%s"       % (self.dir_tophat))
        l_sh_info.append("dir_cufflinks_uk=%s" % (self.dir_cufflinks_unknown))
        l_sh_info.append("""
$samtools_exe view  -H                                                      \\
    $tophat_dir/$samp_name/accepted_hits.genome.sort.bam                    \\
    >$tophat_dir/$samp_name/accepted_hits.header.sam                      &&\\
    
$samtools_exe sort  -n -m 200000000                                         \\
    $tophat_dir/$samp_name/accepted_hits.genome.sort.bam                    \\
    $tophat_dir/$samp_name/accepted_hits.genome.sort_name                 &&\\
    
$samtools_exe view                                                          \\
    -o $tophat_dir/$samp_name/accepted_hits.genome.sort_name.sam            \\
    $tophat_dir/$samp_name/accepted_hits.genome.sort_name.bam             &&\\

[ ! -d $HTS_u_dir/$samp_name ] && mkdir -p $HTS_u_dir/$samp_name

$py_exe $deseq_exe                                                          \\
    -s no -f sam -a 10                                                      \\
    $tophat_dir/$samp_name/accepted_hits.genome.sort_name.sam  $unknown_GTF \\
    >$HTS_u_dir/$samp_name/$samp_name.dexseq_NeoRaw.txt
        """)
        return l_sh_info
    
    def s_05_1_mergeGTF(self):
        l_sh_info    = []
        
        dir_db = "%s/%s" % (self.Database, self.ref)
        self.merge_GTF      = "%s/all.exon.sort.gtf"      % (dir_db)
        self.merge_GTF_ercc = "%s/all.exon.sort.ERCC.gtf" % (dir_db)
        
        l_sh_info.append("unknown_GTF=$1")
        l_sh_info.append("extra_GTF=$2")
        l_sh_info.append("py_exe=%s"        % (self.sftw_py))
        l_sh_info.append("known_GTF=%s"     % (self.genome_gtf))
        l_sh_info.append("merge_GTF=%s"     % (self.merge_GTF))
        l_sh_info.append("merge_ERCC_GTF=%s"% (self.merge_GTF_ercc))
        l_sh_info.append("AddExtra=%s/s05_addExtraInfo.py"  % (self.bin))
        l_sh_info.append("""
sed 's/XLOC_/novoXLOC_/g'        $unknown_GTF                              |\\
     sed 's/TCONS_/novoTCONS_/g' /dev/stdin   >$unknown_GTF.tmp           &&\\

if [ $extra_GTF == "None" ]
    then grep -P "^chr"               $known_GTF                           |\\
             cat /dev/stdin           $unknown_GTF.tmp                     |\\
             bedtools sort -i /dev/stdin                                   |\\
             grep -P "^chr" >$merge_GTF                                   &&\\
         rm $unknown_GTF.tmp                                              &&\\
         grep -P "^ERCC|^RGC" $known_GTF                                   |\\
         cat $merge_GTF /dev/stdin >$merge_ERCC_GTF
fi

if [ $extra_GTF != "None" ]
    then $py_exe $AddExtra $extra_GTF ${extra_GTF}.tmp lncGene lncTran    &&\\
         grep -P "^chr"               $known_GTF                           |\\
             cat /dev/stdin $unknown_GTF.tmp  ${extra_GTF}.tmp             |\\
             bedtools sort -i /dev/stdin                                   |\\
             grep -P "^chr" >$merge_GTF                                   &&\\
         rm $unknown_GTF.tmp ${extra_GTF}.tmp                             &&\\
         grep -P "^ERCC|^RGC" $known_GTF                                   |\\
         cat $merge_GTF /dev/stdin >$merge_ERCC_GTF
fi
        """)
        return l_sh_info
        

    def s_05_2_erccGTF(self, given_GTF, given_GTF_ERCC):
        l_sh_info    = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("given_GTF=%s"        % (given_GTF))
        l_sh_info.append("given_GTF_ERCC=%s"   % (given_GTF_ERCC))
        l_sh_info.append("dir_database=%s/$ref"% (self.Database))
        l_sh_info.append("dir_path=%s"         % (self.path))
        l_sh_info.append("""
cat $given_GTF $dir_path/database/ERCC.gtf >$dir_database/$given_GTF_ERCC
        """)
        return l_sh_info

        
    def s_06_cfqt(self, cflk_dir):
        l_sh_info = []
        l_sh_info.append("brief_name=$1")
        l_sh_info.append("use_gtf=$2")
        l_sh_info.append("sftw_cflk_dir=%s" % (self.sftw_cflk_dir))
        l_sh_info.append("tophat_dir=%s"    % (self.dir_tophat))
        l_sh_info.append("dir_cufflinks=%s" % (cflk_dir))
        l_sh_info.append("""
in_bam=$tophat_dir/$brief_name/accepted_hits.bam
out_dir=$dir_cufflinks/$brief_name

$sftw_cflk_dir/cuffquant                                                    \\
   -p 8  -u                                                                 \\
   -o $out_dir                                                              \\
   $use_gtf                                                                 \\
   $in_bam
        """)
        return l_sh_info
        
    def s_07_cfnm(self, cfnm_dir):
        l_sh_info = []
        l_sh_info.append("list_brief=$1")
        l_sh_info.append("use_gtf=$2")
        l_sh_info.append("sftw_cflk_dir=%s" % (self.sftw_cflk_dir))
        l_sh_info.append("dir_cuffnorm=%s"  % (cfnm_dir))
        l_sh_info.append("""
shift
shift
$sftw_cflk_dir/cuffnorm                                                     \\
    -p 8  -o $dir_cuffnorm  -L $list_brief  $use_gtf                        \\
    $@
        """)
        return l_sh_info


    def s_06_1_cflk(self, cflk_dir):
        l_sh_info = []
        l_sh_info.append("brief_name=$1")
        l_sh_info.append("use_gtf=$2")
        l_sh_info.append("sftw_cflk_dir=%s" % (self.sftw_cflk_dir))
        l_sh_info.append("tophat_dir=%s"    % (self.dir_tophat))
        l_sh_info.append("dir_cufflinks=%s" % (cflk_dir))
        l_sh_info.append("py_exe=%s"           % (self.sftw_py))
        l_sh_info.append("py_grep_FPKM=%s/grep_FPKM.py"  % (self.bin))
        l_sh_info.append("""
in_bam=$tophat_dir/$brief_name/accepted_hits.bam
out_dir=$dir_cufflinks/$brief_name

#$sftw_cflk_dir/cufflinks                                                    \\
#   -p 8                                                                     \\
#   -o $out_dir                                                              \\
#   -G $use_gtf                                                              \\
#   $in_bam

$py_exe $py_grep_FPKM $out_dir/genes.fpkm_tracking $use_gtf                 \\
    >$out_dir/$brief_name.FPKM.txt
        """)
        return l_sh_info


    def s_08_repeatCount(self):
        l_sh_info    = []
        l_sh_info.append("brief_name=$1")
        l_sh_info.append("samtools_exe=%s" % (self.sftw_samtools))
        l_sh_info.append("bedtools_exe=%s" % (self.sftw_bedtools))
        l_sh_info.append("tophat_dir=%s"   % (self.dir_tophat))
        l_sh_info.append("repeat_dir=%s"   % (self.dir_repeat_counts))
        l_sh_info.append("rep_bed=%s"      % (self.rmsk_bed))
        l_sh_info.append("py_exe=%s"       % (self.sftw_py))
        l_sh_info.append("cnt_py=%s/Repeat_Intersect2Count.py"  % (self.bin) )
        
        
        l_sh_info.append("dir_cufflinks_uk=%s" % (self.dir_cufflinks_unknown))
        l_sh_info.append("""
in_bam=$tophat_dir/$brief_name/accepted_hits.bam
out_dir=$repeat_dir/$brief_name

$samtools_exe view -F 0x0004 $in_bam                                       |\\
    grep -v ERCC-00* | grep -v RGC-CRE                                     |\\
    grep -v RGC-GFP  | grep -v RGC-mRFP |grep '\bNH:i:1\b'                 |\\
    awk '{OFS="\\t"; print $3,$4,$4+length($10),$1 }'                       \\
       >${out_dir}/repeat_result.bed

$bedtools_exe intersect -sorted -loj                                        \\
    -a $rep_bed  -b ${out_dir}/repeat_result.bed                           |\\
    $py_exe $cnt_py /dev/stdin >${out_dir}/repeat_count.bed
        """)
        return l_sh_info
        
    def stat_02_SplitHTS(self):
        l_sh_info    = []
        
        dir_HTS = self.dir_HTS_result
        infile     = "%s/merge.dexseq_clean.gene.xls"        % (dir_HTS)
        out_Refseq = "%s/merge.dexseq_clean_refseq.gene.xls" % (dir_HTS)
        out_lncRNA = "%s/merge.dexseq_clean_lncRNA.gene.xls" % (dir_HTS)
        inNeo      = "%s/merge.dexseq_NeoRaw.gene.xls"       % (dir_HTS)
        outNeo     = "%s/merge.dexseq_NeoPass.gene.xls"      % (dir_HTS)
        
        l_sh_info.append("NeoINFO=$1")
        l_sh_info.append("infile=%s"     % (infile    ))
        l_sh_info.append("out_Refseq=%s" % (out_Refseq))
        l_sh_info.append("out_lncRNA=%s" % (out_lncRNA))
        l_sh_info.append("inNeo=%s"      % (inNeo     ))
        l_sh_info.append("outNeo=%s"     % (outNeo    ))
        l_sh_info.append("""
grep -v -P '^lncGene|novoXLOC_' $infile >$out_Refseq
head -n 1 $infile >$out_lncRNA && grep -P '^lncGene' $infile >>$out_lncRNA

head -n 1 $inNeo >$outNeo
for i in `cut -f 1 $NeoINFO | uniq`;do grep -w $i $inNeo ;done >>$outNeo
        """)
        return l_sh_info
        
