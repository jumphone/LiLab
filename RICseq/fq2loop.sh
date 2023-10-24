R1=/home/database/data/FIST2021/data/10k_SRR8179797_1.fastq
R2=/home/database/data/FIST2021/data/10k_SRR8179797_2.fastq
OUT=/home/database/RICseq_demo/test2
CPU=20
RESOLUTION=1000

########
# Human
REF=/home/database/reference/hg19/hg19
CHR_SIZE=/home/database/annotation/hg19/hg19_pure.chrom.sizes
########

##############################################################
# Staring
mkdir $OUT

##############################################################
# QC & Mapping ...
trim_galore --paired $R1 $R2 --length 20 -o $OUT --basename trim_galore_output
VR1=$OUT\/trim_galore_output_val_1.fq
VR2=$OUT\/trim_galore_output_val_2.fq
bowtie2 -x $REF --threads $CPU -U $VR1 --reorder | samtools view -Shb - > $OUT\/read1.bam
bowtie2 -x $REF --threads $CPU -U $VR2 --reorder | samtools view -Shb - > $OUT\/read2.bam
samtools sort $OUT\/read1.bam -o $OUT\/read1.sorted.bam
samtools sort $OUT\/read2.bam -o $OUT\/read2.sorted.bam
samtools rmdup -s $OUT\/read1.sorted.bam $OUT\/read1.sorted.rmdup.bam
samtools rmdup -s $OUT\/read2.sorted.bam $OUT\/read2.sorted.rmdup.bam
samtools merge $OUT\/merged.bam $OUT\/read1.sorted.rmdup.bam $OUT\/read2.sorted.rmdup.bam
samtools sort -n $OUT\/merged.bam -o $OUT\/merged.paired.bam --threads $CPU
bedtools bamtobed -bedpe -i $OUT\/merged.paired.bam >  $OUT\/merged.paired.bedpe

##############################################################
# Formatting ...
cooler cload pairs -c1 1 -p1 2 -c2 4 -p2 5 $CHR_SIZE\:$RESOLUTION   $OUT\/merged.paired.bedpe  $OUT\/mat_$RESOLUTION\.cool
hicConvertFormat --matrices  $OUT\/mat_$RESOLUTION\.cool  --outFileName  $OUT\/mat_$RESOLUTION\.h5 \
                 --inputFormat cool --outputFormat h5
hicConvertFormat --matrices  $OUT\/mat_$RESOLUTION\.h5  --outFileName  $OUT\/mat_$RESOLUTION\.hicpro.mat \
                 --inputFormat h5 --outputFormat hicpro --bedFileHicpro $OUT\/mat_$RESOLUTION\.hicpro.bed

##############################################################
# Loop calling ...
mkdir $OUT\/mat_$RESOLUTION\.fithic
/home/toolkit/local/bin/hicpro2fithic.py  -i $OUT\/mat_$RESOLUTION\.hicpro.mat -b $OUT\/mat_$RESOLUTION\.hicpro.bed -o $OUT\/mat_$RESOLUTION\.fithic -r $RESOLUTION
/home/toolkit/local/bin/fithic -i $OUT\/mat_$RESOLUTION\.fithic/fithic.interactionCounts.gz -f $OUT\/mat_$RESOLUTION\.fithic/fithic.fragmentMappability.gz -o $OUT\/mat_$RESOLUTION\.fithic -r $RESOLUTION -x All
