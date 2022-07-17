#include "MatchBoxPC.h"

void PROCESS_CROSS_EDGE(vector<MilanLongInt> &Counter,
                        MilanLongInt edge,
                        MilanLongInt *SPtr)
{
    // Start: PARALLEL_PROCESS_CROSS_EDGE_B
    MilanLongInt captureCounter;

#pragma omp atomic capture
    captureCounter = --Counter[edge]; // Decrement

    if (captureCounter == 0)
#pragma omp atomic
        (*SPtr)--; // Decrement S

#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << ")Decrementing S: Ghost vertex " << edge << " has received all its messages";
    fflush(stdout);
#endif

    // End: PARALLEL_PROCESS_CROSS_EDGE_B
}