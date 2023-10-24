


# https://liz-fernandez.github.io/HiC-Langebio/04-matrix.html

# https://github.com/4dn-dcic/pairix/blob/master/util/bam2pairs/README.md

RESOLUTION=1000
CHR_SIZE=/home/database/annotation/hg19/hg19_pure.chrom.sizes
INPUT_BAM=/home/xiaodong/ATAC/demo1/Mapping/ESC_G1.pe.q10.sort.rmdup.bam
OUTPUT_PREFIX=/home/database/repli/RICseq_demo/test1

CPU=10

samtools sort -n $INPUT_BAM -o $OUTPUT_PREFIX\_nsort.bam --threads $CPU
bedtools bamtobed -bedpe -i $OUTPUT_PREFIX\_nsort.bam > $OUTPUT_PREFIX\_nsort.bedpe
cooler cload pairs -c1 1 -p1 2 -c2 4 -p2 5    $CHR_SIZE\:$RESOLUTION   $OUTPUT_PREFIX\_nsort.bedpe  $OUTPUT_PREFIX\_$RESOLUTION\.cool

hicConvertFormat --matrices  $OUTPUT_PREFIX\_$RESOLUTION\.cool  --outFileName  $OUTPUT_PREFIX\_$RESOLUTION\.h5 \
                 --inputFormat cool --outputFormat h5


hicConvertFormat --matrices  $OUTPUT_PREFIX\_$RESOLUTION\.cool  --outFileName  $OUTPUT_PREFIX\_$RESOLUTION\.hicpro.mat \
                 --inputFormat cool --outputFormat hicpro --bedFileHicpro $OUTPUT_PREFIX\_$RESOLUTION\.hicpro.bed


mkdir $OUTPUT_PREFIX\_$RESOLUTION\.fithic
/home/toolkit/local/bin/hicpro2fithic.py  -i $OUTPUT_PREFIX\_$RESOLUTION\.hicpro.mat -b $OUTPUT_PREFIX\_$RESOLUTION\.hicpro.bed -o $OUTPUT_PREFIX\_$RESOLUTION\.fithic -r $RESOLUTION

/home/toolkit/local/bin/fithic -i $OUTPUT_PREFIX\_$RESOLUTION\.fithic/fithic.interactionCounts.gz -f $OUTPUT_PREFIX\_$RESOLUTION\.fithic/fithic.fragmentMappability.gz -o $OUTPUT_PREFIX\_$RESOLUTION\.fithic -r $RESOLUTION -x All

