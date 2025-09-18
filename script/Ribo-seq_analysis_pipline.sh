#step1
#trim adaptor
fastp -w 15 -i Mock-48h-1_1.fq.gz -I Mock-48h-1_2.fq.gz  -m --merged_out ./$output"_merge_clean.fastq.gz" --json $output".json" --html $output".html" -q 20 -u 40 -n 0 \
--length_required 23 --length_limit 36 --detect_adapter_for_pe --correction --overlap_len_require 23 --overlap_diff_limit 2 --include_unmerged -r
done


#step2
#remove rRNA/tRNA, download from:https://rnacentral.org/
#BSUB -n 12
#BSUB -o rRNA.log
source ~/env/samtools.env
ls ../*gz |while read id 
do
bowtie2 -x /gss1/home/zpchen/ref/rice/rRNA3/all_remove_rna -U $id --un-gz $(basename $id .fastq.gz)"_rRNA.fastq.gz" -p 12 |samtools view \
-@ 6 -bS > $(basename $id .fastq.gz)".bam"
done


#step3
#align
#BSUB -n 12
#BSUB -o star.log
ls ../*.gz |while read id
do 
STAR --outFilterType BySJout --runThreadN 12 --outFilterMismatchNmax 2 --genomeDir /gss1/home/zpchen/ref/rice/ribo_star \
--readFilesCommand zcat --readFilesIn $id --outFileNamePrefix $(basename $id _merge_clean_rRNA.fastq.gz) --outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 1 --alignEndsType EndToEnd --outSAMattributes All
done

#step4
#Assess ribosome occupancy of lncRNA-derived sORFs
#BSUB -o count.log
cat config |while read id
do
bedtools coverage -split -a predict_rm_overlap_sORF.bed -b $id >$(basename $id .bam)".txt"
done


