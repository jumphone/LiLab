# Reference genome must be hg19 !!!

############################################################
# step01. copy PIQ folder to your path
OLD_PATH=/home/database/download/thashim-piq-single-0e1036763bcc/
PIQ_PATH=/home/database/tmp/PIQ

cp -r $OLD_PATH $PIQ_PATH


####################################
# step02. generate motif Rdata

PIQ_PATH=/home/database/tmp/PIQ
RSCR=/home/database/download/R-3.1.2/bin/Rscript

cd $PIQ_PATH; $RSCR pwmmatch.exact.r common.r  pwms/jaspar_select.txt 1 ./motif.matches/
cd $PIQ_PATH; $RSCR pwmmatch.exact.r common.r  pwms/jaspar_select.txt 2 ./motif.matches/


###################################
# step03. use PIQ to pre-process your BAM file

PIQ_PATH=/home/database/tmp/PIQ
RSCR=/home/database/download/R-3.1.2/bin/Rscript

BAM=/home/database/tmp/tmp.bam
OUTPUT=/home/database/tmp/tmp.bam.PIQ

mkdir $OUTPUT

cd $PIQ_PATH; $RSCR pairedbam2rdata.r common.r $OUTPUT\/BAM.RData $BAM 


#####################################
# step04. run PIQ

PIQ_PATH=/home/database/tmp/PIQ
RSCR=/home/database/download/R-3.1.2/bin/Rscript

BAM=/home/database/tmp/tmp.bam
OUTPUT=/home/database/tmp/tmp.bam.PIQ

TAG=TFAP2C
NUM=1
###########
TMP=$OUTPUT\/$TAG\.PIQ.TMP
OUT=$OUTPUT\/$TAG\.PIQ.OUT
mkdir $OUT
cd $PIQ_PATH; $RSCR pertf.r common.r ./motif.matches/  $TMP  $OUT  $OUTPUT\/BAM.RData  $NUM






