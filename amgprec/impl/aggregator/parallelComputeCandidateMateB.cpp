#include "MatchBoxPC.h"

void PARALLEL_COMPUTE_CANDIDATE_MATE_B(MilanLongInt NLVer,
                                              MilanLongInt *verLocPtr,
                                              MilanLongInt *verLocInd,
                                              MilanInt myRank,
                                              MilanReal *edgeLocWeight,
                                              MilanLongInt *candidateMate)
{

    MilanLongInt v = -1;

#pragma omp parallel private(v) default(shared) num_threads(NUM_THREAD)
    {

#pragma omp for schedule(static)
        for (v = 0; v < NLVer; v++) {
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
