rm amgprec/impl/aggregator/algoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMateCMP.o
make all
cd samples/advanced/pdegen
make amg_d_pde3d
cd runs
mpirun -np 4 amg_d_pde3d amg_pde3d.inp



