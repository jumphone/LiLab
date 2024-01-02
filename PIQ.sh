
PIQ_PATH=/home/database/download/thashim-piq-single-0e1036763bcc/
RSCR=/home/database/download/R-3.1.2/bin/Rscript

cd $PIQ_PATH; $RSCR pwmmatch.exact.r common.r  pwms/jaspar_select.txt 1 ./motif.matches/
cd $PIQ_PATH; $RSCR pwmmatch.exact.r common.r  pwms/jaspar_select.txt 2 ./motif.matches/
cd $PIQ_PATH; $RSCR pwmmatch.exact.r common.r  pwms/jaspar_select.txt 3 ./motif.matches/
cd $PIQ_PATH; $RSCR pwmmatch.exact.r common.r  pwms/jaspar_select.txt 4 ./motif.matches/


#cd $PIQ_PATH; $RSCR pwmmatch.exact.r common.r  pwms/jaspar_select.txt 5 ./motif.matches/


:<<!
BAM=/home/database/data/ATAC_CLL/X101SC21013172-Z01-J005/1.rawdata/OUTPUT_LOG/TROP2_Hi.bam
cd $PIQ_PATH; nohup $RSCR pairedbam2rdata.r common.r $BAM\.PIQ.RData $BAM &

BAM=/home/database/data/ATAC_CLL/X101SC21013172-Z01-J005/1.rawdata/OUTPUT_LOG/TROP2.bam
cd $PIQ_PATH; nohup $RSCR pairedbam2rdata.r common.r $BAM\.PIQ.RData $BAM &
!


