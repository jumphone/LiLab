

# https://github.com/4dn-dcic/pairix/blob/master/util/bam2pairs/README.md

bam2pairs -l -c CHR_SIZE.txt    INPUT.bam    OUTPUT_PREFIX

cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5    CHR_SIZE.txt:10000     OUTPUT_PREFIX.bsorted.pairs     OUTPUT_PREFIX_10k.cool


