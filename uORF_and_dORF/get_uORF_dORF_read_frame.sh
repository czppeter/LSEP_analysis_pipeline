
#step1: get 5UTR and 3UTR fasta sequence
cat rice_5UTR.txt |awk -F"\t" 'BEGIN{FS=OFS="\t"}{print $1,$4-1,$5,$3,$7,$9}' |cut -d';' -f 1 |awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$6,$4,$5}' >utr_5.bed
bedtools getfasta -fi ../tigr7.fa -bed utr_5.bed -s -name >utr_5_2.fa

cat rice_3UTR.txt |awk -F"\t" 'BEGIN{FS=OFS="\t"}{print $1,$4-1,$5,$3,$7,$9}' |cut -d';' -f 1 |awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$6,$4,$5}' >utr_3.bed
bedtools getfasta -fi ../tigr7.fa -bed utr_3.bed -s -name >utr_3.fa

less utr_5.fa |grep '>' |sed 's/:/\t/' |sed 's/://' |awk '{a[$1]=a[$1]";"$2}END{for(i in a){print i"\t"a[i]}}' |sort -k1,1 |sed 's/;//'  >utr_5_pos_message.txt
less -SN utr_5.fa |grep '>' |cut -d':' -f 1 |cut -d'=' -f 2 >1
cat utr_5.fa |grep -v '>' >2
paste 1 2 |awk '{seq[$1]=seq[$1]$2}END{for(i in seq ){print ">"i";utr5";print seq[i]}}' >utr_5_all_seq.fa


#step2: get uORF sequence
conda activate emboss
getorf -sequence utr_5_all_seq.fa -table 1 -minsize 15 -maxsize 450 -find 1 -outseq uORF_seq/uORF_450.fa
getorf -sequence utr_5_all_seq.fa -table 1 -minsize 15 -maxsize 450 -find 3 -outseq uORF_seq/uORF_450_3.fa

awk 'NR>1 && /^>/ {print ""}{printf "%s", $0":"}' uORF_450_3.fa |cut -d':' -f 2- |sed 's/://g'  >1
grep '>' uORF_450_3.fa >2
paste 2 1 |grep -v 'REVERSE SENSE' >uORF_450_nt.fa

# get uORF end codon position
less -SN uORF_450_3.fa |grep '>' |grep -v 'REVERSE SENSE' |cut -d'-' -f 2 |sed 's/]//g'  >1
less -SN uORF_450_3.fa |grep '>' |grep -v 'REVERSE SENSE' |cut -d'_' -f 1,2  >2
paste 2 1 >end_codon_pos.txt

# get uORF sequence with end codon
awk 'ARGIND==1{a[$1]=$2}ARGIND>1 {print $0"\t" substr(a[$1],$2+1,3)}' ../utr_5_all_seq.txt end_codon_pos.txt >end_codon.txt
paste uORF_450_nt.fa end_codon.txt |cut -f 1,2,5 |awk -F"\t" '{if($3=="TAA"||$3=="TAG"||$3=="TGA")print $0}' >all_uORF_seq.txt

# get the genomic regions of uORFs
cat utr_5_pos_message.txt |sed -E 's/utr_[0-9]+:Chr[0-9]+://g' >utr_5_pos.txt
less -SN all_uORF_seq.txt |cut -f 1 |sed 's/ //g' |sed 's/;/\t/g' |sed 's/\[/\t/g' |sed 's/\]//g' |cut -f 1,3 |sed 's/>//g' >all_uORF_region.txt

cat utr_5_pos.txt |sed -E 's/\([^)]*\)//g'  >1
cat utr_5_pos.txt |cut -d';' -f 1 |cut -d'(' -f 2 |sed 's/)//g' >2
paste 1 2 >utr_5_pos2.txt

#run python get_uORF_bed12.py 


#step3: get 3'UTR sequence
less utr_3.fa |grep '>' |sed 's/:/\t/' |sed 's/://' |awk '{a[$1]=a[$1]";"$2}END{for(i in a){print i"\t"a[i]}}' |sort -k1,1 |sed 's/;//'  >utr_3_pos_message.txt
less -SN utr_3.fa |grep '>' |cut -d':' -f 1 |cut -d'=' -f 2 >1
cat utr_3.fa |grep -v '>' >2
paste 1 2 |awk '{seq[$1]=seq[$1]$2}END{for(i in seq ){print ">"i";utr3";print seq[i]}}' >utr_3_all_seq.fa

#step4: get dORF sequence
conda activate emboss
getorf -sequence utr_3_all_seq.fa -table 1 -minsize 15 -maxsize 450 -find 1 -outseq dORF_seq/dORF_450.fa
getorf -sequence utr_3_all_seq.fa -table 1 -minsize 15 -maxsize 450 -find 3 -outseq dORF_seq/dORF_450_3.fa

awk 'NR>1 && /^>/ {print ""}{printf "%s", $0":"}' dORF_450_3.fa |cut -d':' -f 2- |sed 's/://g'  >1
grep '>' dORF_450_3.fa >2
paste 2 1 |grep -v 'REVERSE SENSE' >dORF_450_nt.fa

# get dORF end codon position
less -SN dORF_450_3.fa |grep '>' |grep -v 'REVERSE SENSE' |cut -d'-' -f 2 |sed 's/]//g'  >1
less -SN dORF_450_3.fa |grep '>' |grep -v 'REVERSE SENSE' |cut -d'_' -f 1,2  >2
paste 2 1 >end_codon_pos.txt

# get dORF sequence with end codon
cat utr_3_all_seq.fa |grep '>' >1
cat utr_3_all_seq.fa |grep -v '>' >2
paste 1 2 >utr_3_all_seq.txt

awk 'ARGIND==1{a[$1]=$2}ARGIND>1 {print $0"\t" substr(a[$1],$2+1,3)}' ../utr_5_all_seq.txt end_codon_pos.txt >end_codon.txt
paste dORF_450_nt.fa end_codon.txt |cut -f 1,2,5 |awk -F"\t" '{if($3=="TAA"||$3=="TAG"||$3=="TGA")print $0}' >all_dORF_seq.txt

# get the genomic regions of dORFs
cat utr_3_pos_message.txt |sed -E 's/utr_[0-9]+:Chr[0-9]+://g' >utr_3_pos.txt
less  all_dORF_seq.txt |cut -f 1 |sed 's/ //g' |sed 's/;/\t/g' |sed 's/\[/\t/g' |sed 's/\]//g' |cut -f 1,3 |sed 's/>//g' >all_dORF_region.txt
cat all_dORF_region.txt |awk -F"-" '{print $1"-"$2+3}' >all_dORF_region2.txt 

cat utr_3_pos.txt |sed -E 's/\([^)]*\)//g'  >1
cat utr_3_pos.txt |cut -d';' -f 1 |cut -d'(' -f 2 |sed 's/)//g' >2
paste 1 2 >utr_3_pos2.txt

#run python get_dORF_bed12.py 
