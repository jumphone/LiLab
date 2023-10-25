

.norm1 <- function(x){
    y=(x-min(x))/(max(x)-min(x))
    y[which(is.na(y))]=0
    return(y)
    }





.stepClst <- function(MAT){
    MAT=MAT
    ######################
    print('scaling...')
    SMAT=t(apply(MAT,1,scale))
    rownames(SMAT)=rownames(MAT)
    colnames(SMAT)=colnames(MAT)
    SMAT[which(is.na(SMAT))]=0
    ###################
    print('scoring...')
    TMP=SMAT
    TMP_ROW_RANK=t(apply(-TMP,1,rank,ties.method='random'))
    TMP_FLAG=TMP_ROW_RANK
    TMP_FLAG[which(TMP_ROW_RANK>1)]=0
    ####################
    COL_INDEX=c(1:ncol(TMP))
    ####################
    DDD=COL_INDEX - median(COL_INDEX)
    WWW=  10 ** abs(DDD) * sign(DDD)
    ####################
    TMP_POS=TMP
    TMP_POS[which(TMP<0)]=0
    SCORE1 = .norm1(TMP_POS %*% WWW )*0.99
    ####################
    SCORE2 = TMP_FLAG %*% COL_INDEX
    SCORE = SCORE1 + SCORE2
    SPLIT = as.integer(SCORE)
    #####################
    print('finished!')
    #####################
    OUT=list()
    OUT$smat=SMAT
    OUT$score=SCORE
    OUT$split=SPLIT
    return(OUT)
    }
