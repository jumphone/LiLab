import sys
fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')

l1=fi.readline()
l2=fi.readline()

while l1 !='':
    seq1=l1.rstrip().split('\t')   
    seq2=l2.rstrip().split('\t')   
    this_out=[seq1[0],seq1[1],seq1[2],seq2[0],seq2[1],seq2[2],seq1[3],seq1[4],seq1[5],seq2[3],seq2[4],seq2[5]]
    fo.write('\t'.join(this_out)+'\n')
    l1=fi.readline()
    l2=fi.readline()



fo.close()
fi.close()
