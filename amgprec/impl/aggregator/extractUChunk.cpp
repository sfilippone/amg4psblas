#include "MatchBoxPC.h"

void extractUChunk(
    vector<MilanLongInt> &UChunkBeingProcessed,
    vector<MilanLongInt> &U,
    vector<MilanLongInt> &privateU)
{

    UChunkBeingProcessed.clear();
#pragma omp critical(U)
    {

        if (U.empty() && !privateU.empty()) // If U is empty but there are nodes in private U
        {
            while (!privateU.empty())
                UChunkBeingProcessed.push_back(privateU.back());
            privateU.pop_back();
        }
        else
        {
            for (int i = 0; i < UCHUNK; i++)
            { // Pop the new nodes
                if (U.empty())
                    break;
                UChunkBeingProcessed.push_back(U.back());
                U.pop_back();
            }
        }

    } // End of critical U // End of critical U
}