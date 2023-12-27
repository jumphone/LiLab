SAM=./collect_pair_tag.interaction.sam

samtools view -b $SAM -o $SAM\.bam
bedtools bamtobed -i $SAM\.bam > $SAM\.bam.bed
python /home/database/pipeline/RIC/RICpipe/bedtobedpe.py $SAM\.bam.bed $SAM\.bam.bedpe

