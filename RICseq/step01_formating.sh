


# https://liz-fernandez.github.io/HiC-Langebio/04-matrix.html

# https://github.com/4dn-dcic/pairix/blob/master/util/bam2pairs/README.md

RESOLUTION=10000
CHR_SIZE=
INPUT_BAM=
OUTPUT_PREFIX=

bam2pairs -l -c $CHR_SIZE    $INPUT_BAM   $OUTPUT_PREFIX
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5    $CHR_SIZE\:$RESOLUTION   $OUTPUT_PREFIX\.bsorted.pairs  $OUTPUT_PREFIX\_10k.cool


