#/home/toolkit/tools/piswap


from __future__ import print_function
import PSB2009_3opt_2 as psb09
import matching as match
G = psb09.getGraph("network1.tab") #input network1
G2 = psb09.getGraph("network2.tab")#input network2
GS = psb09.graphScores("pairwise_sequence_similarity_of_network1_and_2.evals")


M0 = match.max_weight_matching(GS) #run hungorian algorith
(S, M) = psb09.processOnce(G, G2, GS, M0, 0.6, 200) #run PISwap

F = open("match_output.txt","w")
for node in M:
    F.write(node+" "+M[node]+'\n');

F.close()
