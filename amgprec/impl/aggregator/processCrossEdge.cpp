#include "MatchBoxPC.h"

void PROCESS_CROSS_EDGE(vector<MilanLongInt> &Counter,
                               map<MilanLongInt, MilanLongInt> &Ghost2LocalMap,
                               MilanLongInt edge,
                               MilanLongInt *SPtr)
{
    MilanLongInt S = *SPtr;
    // Decrement the counter:
    // Start: PARALLEL_PROCESS_CROSS_EDGE_B
    if (Counter[Ghost2LocalMap[edge]] > 0)
    {
        Counter[Ghost2LocalMap[edge]] -= 1; // Decrement
        if (Counter[Ghost2LocalMap[edge]] == 0)
        {
            S--; // Decrement S
#ifdef PRINT_DEBUG_INFO_
            cout << "\n(" << myRank << ")Decrementing S: Ghost vertex " << edge << " has received all its messages";
            fflush(stdout);
#endif
        }
    } // End of if Counter[edge] > 0
      // End: PARALLEL_PROCESS_CROSS_EDGE_B
    *SPtr = S;
}