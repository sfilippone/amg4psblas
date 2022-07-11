#include "MatchBoxPC.h"

/**
 * //TODO documentation
 * @param k
 * @param verLocInd
 * @param StartIndex
 * @param EndIndex
 * @param GMate
 * @param Mate
 * @param Ghost2LocalMap
 * @return
 */
bool isAlreadyMatched(MilanLongInt node,
                             MilanLongInt StartIndex,
                             MilanLongInt EndIndex,
                             vector <MilanLongInt> &GMate,
                             MilanLongInt* Mate,
                             map <MilanLongInt, MilanLongInt> &Ghost2LocalMap
) {

    bool result = false;
#pragma omp critical(Mate)
    {
        if ((node < StartIndex) || (node > EndIndex)) { //Is it a ghost vertex?
            if (GMate[Ghost2LocalMap[node]] >= 0)// Already matched
                result = true;
        } else { //A local vertex
            if (Mate[node - StartIndex] >= 0) // Already matched
                result = true;
        }

    }

    return result;
}