cd amgprec/impl/aggregator/
rm MatchBoxPC.o 
rm sendBundledMessages.o 
rm initialize.o 
rm extractUChunk.o 
rm isAlreadyMatched.o 
rm findOwnerOfGhost.o 
rm computeCandidateMate.o 
rm parallelComputeCandidateMateB.o 
rm processMatchedVertices.o 
rm processCrossEdge.o 
rm queueTransfer.o 
rm processMessages.o 
rm processExposedVertex.o 
rm algoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMateC.o 
rm algoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMateCMP.o
cd ../../../
make all
cd samples/advanced/pdegen
make amg_d_pde3d
cd runs
mpirun -np 4 amg_d_pde3d amg_pde3d.inp



