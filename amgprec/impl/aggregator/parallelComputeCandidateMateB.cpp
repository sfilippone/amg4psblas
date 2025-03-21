#include "amg_config.h"
#include "MatchBoxPC.h"
#if !defined(PSB_SERIAL_MPI)

void PARALLEL_COMPUTE_CANDIDATE_MATE_BD(MilanLongInt NLVer,
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
            candidateMate[v] = firstComputeCandidateMateD(verLocPtr[v], verLocPtr[v + 1],
							 verLocInd, edgeLocWeight);
            // End: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
        }
    }
}

void PARALLEL_COMPUTE_CANDIDATE_MATE_BS(MilanLongInt NLVer,
                                              MilanLongInt *verLocPtr,
                                              MilanLongInt *verLocInd,
                                              MilanInt myRank,
                                              MilanFloat *edgeLocWeight,
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
            candidateMate[v] = firstComputeCandidateMateS(verLocPtr[v], verLocPtr[v + 1],
							 verLocInd, edgeLocWeight);
            // End: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
        }
    }
}
#endif
