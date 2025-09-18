#step1
#trim adaptor
trim_galore -j 15 -q 25 --phred33 --length 35 -e 0.1 --stringency 2 --paired sample_R1.fastq.gz  sample_R2.fastq.gz -o clean_data/

#step2
#align fastq file 
#config file contain sample file name
#for example:Mock-48h-1-lnc_R1_val_1.fq.gz	Mock-48h-1-lnc_R2_val_2.fq.gz

#BSUB -n 15
#BSUB -o align2.log
source ~/env/samtools.env
cat ../config |while read id
do
arr=($id)
fq1=${arr[0]}
fq2=${arr[1]}
output=${fq1%%_*}
hisat2 --dta -x ~/ref/rice/hisat2_rna_seq/tigr7 --rna-strandness RF -p 15  -1 ../$fq1 -2 ../$fq2 --new-summary --summary-file $output'.ht2.txt' |samtools view -@ 15 -bS  |samtools sort -@ 15 -o $output'.sort.bam'
samtools index $output'.sort.bam'
samtools flagstat  $output'.sort.bam' > $output'.sort.stat'

samtools view -@ 15 -q 30 -o  $output'.sort.q30.bam'  $output'.sort.bam'
samtools index  $output'.sort.q30.bam'
samtools flagstat $output'.sort.q30.bam' > $output'.sort.q30.stat'
done

#step3
#config file contain sample bam file;for example:Mock-48h-1-lnc.uniq.bam
#BSUB -o stringtie.log
#BSUB -n 8
cat ../config |while read id
do
output=${id%%.*}
stringtie -p 8 -G /gss1/home/zpchen/ref/rice/tigr_gff/all.gff3 -o $output".gtf" -l $output ../$id
done

stringtie --merge -p 8 -G /gss1/home/zpchen/ref/rice/tigr_gff/all.gff3 -o all_merged.gtf -l lncrna gtf.list.txt
gffcompare -r /gss1/home/zpchen/ref/rice/tigr_gff/all.gff3 -o all_merged_compare.gtf all_merged.gtf

#get "i、u、x、o" transcripit
awk '$3=="i" || $3=="u" || $3=="x" || $3=="o" || $3=="j" && $10>200 {print "\"" $5 "\""}' all_merged_compare.gtf.all_merged.gtf.tmap >filter_iuxoj_trans.txt
fgrep -f filter_iuxoj_trans.txt all_merged_compare.gtf.annotated.gtf >filter_iuxoj_trans.gtf

# get transcript sequence
gffread filter_iuxoj_trans.gtf -g  /gss1/home/zpchen/ref/rice/tigr7.fa -w filter_iuxoj_trans.fa

#stpe4 Prediction of lncRNA coding potential (CPC2, LGC, Pfam_scan)
conda activate cpc2

#BSUB -o cpc2.log
ls *.fa |while read id
do
CPC2.py -i filter_iuxoj_trans.fa -o filter_iuxoj_trans_cpc2.out --ORF
done

#BSUB -o LGC.log
source ~/env/python2.env
ls *.fa |while read id
do
python /gss1/home/zpchen/biosoft/LGC/LGC-1.0.py filter_iuxoj_trans.fa filter_iuxoj_trans_LGC.out
done

conda activate emboss

ls *.fa |while read id
do
transeq filter_iuxoj_trans.fa filter_iuxoj_trans.pep -frame=6
done

conda activate pfam_scan
#BSUB -o pfam_scan2.log
#BSUB -n 10
ls *.pep |while read id
do
pfam_scan.pl -cpu 10 -fasta  filter_iuxoj_trans.pep -dir /gss1/home/zpchen/miniconda3/envs/pfam_scan/pfam/ -outfile  filter_iuxoj_trans.pep_pfam.out
done

# step5 retained non-coding lncRNAs that appeared in all prediction tools
cat filter_iuxoj_trans_cpc2.out.txt |awk '{if($9=="noncoding") print $0}' >cpc2
cat filter_iuxoj_trans_LGC.out |awk '{if($5=="Non-coding")print $0}' >LGC
grep -v '^#' filter_iuxoj_trans_pfam.out | grep -v '^\s*$' | awk '($13< 1e-5){print $1}'| awk -F "_" '{print$1}' | sort | uniq >pfam_coding

grep '>' filter_iuxoj_trans.fa |cut -d ' '  -f 1 |sed 's/>//g' >all_lncRNA.txt
awk 'ARGIND==1{a[$1]=$1} ARGIND>1 && !($1 in a){print $0}' pfam_coding all_lncRNA.txt >pfam_noncoding

cat cpc2 |cut -f 1 |sort -k1,1 |uniq  >1
cat LGC |cut -f 1 |sort -k1,1 |uniq  >2
awk 'ARGIND==1{a[$1]}ARGIND>1 && ($1 in a){print $1}' 1 2 >3
awk 'ARGIND==1{a[$1]}ARGIND>1 && ($1 in a){print $1}'  pfam_noncoding 3  |sort -k1,1 >final_nocoding_trans.txt

cat filter_iuxoj_trans.gtf |cut -f 9 |cut -d ';' -f 1 |cut -d '"' -f 2 >1
paste filter_iuxoj_trans.gtf 1 >2
awk -F"\t" 'ARGIND==1{a[$1]=$1} ARGIND>1 && ($10 in a) {print $0}' final_nocoding_trans.txt 2 >3
cat 3 |cut -f 1-9 >final_nocoding_trans.gtf
cat final_nocoding_trans.gtf |grep -v 'ChrM' >final_nocoding_trans_rmChrm.gtf

#step6 predicted lncRNA-derived sORF

gffread final_nocoding_trans_rmChrm.gtf -g /gss1/home/zpchen/ref/rice/tigr7.fa -w final_nocoding_trans_rmChrm.fa
conda activate emboss

# nt sequence
getorf -sequence final_nocoding_trans_rmChrm.fa -table 1 -minsize 15 -maxsize 450 -find 1 -outseq sORF_450/sORF_450.fa
#amino acid sequence
getorf -sequence final_nocoding_trans_rmChrm.fa -table 1 -minsize 15 -maxsize 450 -find 3 -outseq sORF_450/sORF_450_3.fa
grep '>' sORF_450.fa |grep -v 'REVERSE SENSE' >all_sORF.txt

less -SN final_nocoding_trans.gtf |awk '{if($7==".")print $0}' >all_nostrand.gtf
less -SN all_nostrand.gtf  |cut -d';' -f 1 |cut -f 9 |cut -d'"' -f 2 |sort -k1,1 |uniq >no_strand_trans.txt

awk 'NR>1 && /^>/ {print ""}{printf "%s", $0";"}' sORF_450.fa |cut -d';' -f 2- |sed 's/;//g' >1
grep '>' sORF_450.fa >2
paste 2 1  |grep -v '(REVERSE SENSE)' >3
paste 2 1 |grep '(REVERSE SENSE)' >rev_strand_aa.txt
fgrep -f ../no_strand_trans.txt 2 |grep 'REVERSE SENSE' > all_nostrand_rev.txt


cat 3 |cut -d' ' -f 1 |sed 's/>//g' >1
cat 3 |cut -f 2  >2
paste 1 2  >sORF_450_aa.txt


#A small subset of sequences that cannot be distinguished between the sense and antisense strands were also taken into consideration
awk -F"\t" 'ARGIND==1{a[$1]=$0}ARGIND>1 && ($1 in a){print a[$1]}' rev_strand_aa.txt all_nostrand_rev.txt  >all_nostrand_rev_aa.txt
cat all_nostrand_rev_aa.txt |cut -f 1 |cut -d' ' -f 1 |sed 's/>//g' >1
cat all_nostrand_rev_aa.txt  |cut -f 2  >2
paste 1 2  >all_nostrand_rev_aa2.txt

cat sORF_450_aa.txt all_nostrand_rev_aa2.txt >sORF_450_aa_all.txt


#get sORF bed file
python get_sORF_bed.py

# remove repeat sORF
python get_rm_full_overlap.py












