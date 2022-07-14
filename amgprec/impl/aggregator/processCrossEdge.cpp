#include "MatchBoxPC.h"

void PROCESS_CROSS_EDGE(vector<MilanLongInt> &Counter,
                        MilanLongInt edge,
                        MilanLongInt *SPtr)
{
    // Decrement the counter:
    // Start: PARALLEL_PROCESS_CROSS_EDGE_B
    if (Counter[edge] > 0)
    {
        Counter[edge] -= 1; // Decrement
        if (Counter[edge] == 0)
        {
            (*SPtr)--; // Decrement S
#ifdef PRINT_DEBUG_INFO_
            cout << "\n(" << myRank << ")Decrementing S: Ghost vertex " << edge << " has received all its messages";
            fflush(stdout);
#endif
        }

    } // End of if Counter[edge] > 0
      // End: PARALLEL_PROCESS_CROSS_EDGE_B
}