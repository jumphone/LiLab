

R1=/home/database/data/FIST2021/data/10k_SRR8179797_1.fastq
R2=/home/database/data/FIST2021/data/10k_SRR8179797_2.fastq

OUT=/home/database/pipeline/RIC/output



########################################################################
mkdir $OUT

STAR=/home/database/pipeline/RIC/STAR-2.5.2b/bin/Linux_x86_64/STAR
Trimmomatic=/home/database/pipeline/RIC/Trimmomatic-0.36/trimmomatic-0.36.jar
rRNA=/home/database/reference/rRNA_star
GENOME=/home/database/reference/hg19_star
RICpipe=/home/database/pipeline/RIC/RICpipe

GENE_BED=/home/database/pipeline/RIC/data/gencode.v19.annotation.bed
GENE_JUNCTION_BEDPE=/home/database/pipeline/RIC/data/gencode.v19.annotation.bed.junction.bedpe
minimum_fragment_size=150
java -jar $Trimmomatic PE -phred33 -threads 10 $R1 $R2 $OUT\/read1.clean.pair.fq $OUT\/read1.clean.unpair.fq $OUT\/read2.clean.pair.fq $OUT\/read2.clean.unpair.fq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:7:8:true LEADING:25 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:30

perl $RICpipe\/step0.remove_PCR_duplicates/remove_PCR_duplicates.pl $OUT\/read1.clean.pair.fq $OUT\/read2.clean.pair.fq $OUT\/read1.clean.unpair.fq $OUT\/read2.clean.unpair.fq $OUT > $OUT\/run_step0.sh;cd $OUT;sh run_step0.sh

cutadapt -j 10 -b A{100} -b C{100} -b G{100} -b T{100} -n 3 --minimum-length=30 -e 0.1 -o $OUT\/read1.clean.rmDup.rmPoly.fq $OUT\/read1.clean.rmDup.fq
cutadapt -j 10 -b A{100} -b C{100} -b G{100} -b T{100} -n 3 --minimum-length=30 -e 0.1 -o $OUT\/read2.clean.rmDup.rmPoly.fq $OUT\/read2.clean.rmDup.fq

$STAR --runMode alignReads --genomeDir $rRNA --readFilesIn $OUT\/read1.clean.rmDup.rmPoly.fq --outFileNamePrefix $OUT\/read1_torRNA_ \
--outReadsUnmapped Fastx --outFilterMultimapNmax 100 --outSAMattributes All --alignIntronMin 1 --scoreGapNoncan -4 --scoreGapATAC -4 \
--chimSegmentMin 15 --chimJunctionOverhangMin 15 --limitOutSJcollapsed 10000000 --limitIObufferSize 1500000000 --runThreadN 16 \
--alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --alignSJstitchMismatchNmax 5 -1 5 5 

$STAR --runMode alignReads --genomeDir $rRNA --readFilesIn $OUT\/read2.clean.rmDup.rmPoly.fq --outFileNamePrefix $OUT\/read2_torRNA_ \
--outReadsUnmapped Fastx --outFilterMultimapNmax 100 --outSAMattributes All --alignIntronMin 1 --scoreGapNoncan -4 --scoreGapATAC -4 \
--chimSegmentMin 15 --chimJunctionOverhangMin 15 --limitOutSJcollapsed 10000000 --limitIObufferSize 1500000000 --runThreadN 16 \
--alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --alignSJstitchMismatchNmax 5 -1 5 5

$STAR --runMode alignReads --genomeDir $GENOME --readFilesIn $OUT\/read1_torRNA_Unmapped.out.mate1 --outFileNamePrefix $OUT\/read1_toGenome \
--outReadsUnmapped Fastx --outFilterMultimapNmax 100 --outSAMattributes All --alignIntronMin 1 --scoreGapNoncan -4 --scoreGapATAC -4 \
--chimSegmentMin 15 --chimJunctionOverhangMin 15 --limitOutSJcollapsed 10000000 --limitIObufferSize 1500000000 --runThreadN 16 \
--alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --alignSJstitchMismatchNmax 5 -1 5 5

$STAR --runMode alignReads --genomeDir $GENOME --readFilesIn $OUT\/read2_torRNA_Unmapped.out.mate1 --outFileNamePrefix $OUT\/read2_toGenome \
--outReadsUnmapped Fastx --outFilterMultimapNmax 100 --outSAMattributes All --alignIntronMin 1 --scoreGapNoncan -4 --scoreGapATAC -4 \
--chimSegmentMin 15 --chimJunctionOverhangMin 15 --limitOutSJcollapsed 10000000 --limitIObufferSize 1500000000 --runThreadN 16 \
--alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --alignSJstitchMismatchNmax 5 -1 5 5

perl $RICpipe\/step1.collect_pair_tags/collect_pair_tags.pl $OUT\/read1_toGenomeAligned.out.sam $OUT\/read2_toGenomeAligned.out.sam $OUT\/read1_toGenomeChimeric.out.sam $OUT\/read2_toGenomeChimeric.out.sam 10 $OUT\/collect_pair_tag > $OUT\/run_step1.sh;cd $OUT;sh run_step1.sh

perl $RICpipe\/step2.separate_intra_inter_molecular/separate_intra_inter.pl $OUT\/collect_pair_tag.interaction.sam $GENE_BED > $OUT\/run_step2.sh;cd $OUT;sh run_step2.sh
!

perl $RICpipe\/step3.category_intra_reads/category_intra_reads.pl $OUT\/collect_pair_tag.interaction.sam  $minimum_fragment_size $GENE_BED $GENE_JUNCTION_BEDPE > $OUT\/run_step3.sh;cd $OUT;sh run_step3.sh

cp $RICpipe\/lilab_step01_sam2bedpe.sh $OUT
cp $RICpipe\/lilab_step02_statlen.R $OUT

cd $OUT; sh lilab_step01_sam2bedpe.sh
cd $OUT; Rscript lilab_step02_statlen.R














