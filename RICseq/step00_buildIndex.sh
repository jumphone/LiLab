

STAR=/home/database/pipeline/RIC/STAR-2.5.2b/bin/Linux_x86_64/STAR

$STAR \
--runThreadN 20 \
--runMode genomeGenerate \
--genomeSAindexNbases 6 \
--genomeDir /home/database/reference/rRNA_star \
--genomeFastaFiles /home/database/reference/rRNA/human_rRNA_U13369_1.fasta

$STAR \
--runThreadN 20 \
--runMode genomeGenerate \
--genomeSAindexNbases 6 \
--genomeDir /home/database/reference/hg19_star \
--genomeFastaFiles /home/database/reference/hg19/hg19.fa


:<<!
GTF=/home/database/pipeline/RIC/data/gencode.v19.annotation.gtf
BED=/home/database/pipeline/RIC/data/gencode.v19.annotation.bed

perl /home/database/pipeline/RIC/RICpipe/scripts/gtf_to_bed.pl $GTF > $BED\.ori
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' $BED\.ori >  $BED
perl /home/database/pipeline/RIC/RICpipe/scripts/creat_junction_bed.pl $BED\.ori > $BED\.junction.bedpe 
!


