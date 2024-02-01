#include "MatchBoxPC.h"
#ifdef OMP
void PROCESS_CROSS_EDGE(MilanLongInt *edge,
                        MilanLongInt *S)
{
    // Start: PARALLEL_PROCESS_CROSS_EDGE_B
    MilanLongInt captureCounter;

#pragma omp atomic capture
    captureCounter = --(*edge); // Decrement

    //assert(captureCounter >= 0);

    if (captureCounter == 0)
#pragma omp atomic
        (*S)--; // Decrement S

#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << ")Decrementing S: Ghost vertex " << edge << " has received all its messages";
    fflush(stdout);
#endif

    // End: PARALLEL_PROCESS_CROSS_EDGE_B
}
#endif
