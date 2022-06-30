#include "MatchBoxPC.h"
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <map>
#include <vector>
#include "primitiveDataTypeDefinitions.h"
#include "dataStrStaticQueue.h"
#include "omp.h"

inline void PARALLEL_COMPUTE_CANDIDATE_MATE_B(MilanLongInt NLVer,
                                              MilanLongInt *verLocPtr,
                                              MilanLongInt *verLocInd,
                                              MilanInt myRank,
                                              MilanReal *edgeLocWeight,
                                              MilanLongInt *candidateMate)
{

    MilanLongInt v = -1;

#pragma omp parallel private(v) default(shared) num_threads(4)
    {

#pragma omp for schedule(static)
        for (v = 0; v < NLVer; v++)
        {
#ifdef PRINT_DEBUG_INFO_
            cout << "\n(" << myRank << ")Processing: " << v + StartIndex << endl;
            fflush(stdout);
#endif
            // Start: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
            candidateMate[v] = firstComputeCandidateMate(verLocPtr[v], verLocPtr[v + 1], verLocInd, edgeLocWeight);
            // End: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
        }
    }
}
