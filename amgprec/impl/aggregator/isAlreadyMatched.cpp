#include "MatchBoxPC.h"
#ifdef OMP
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
                      vector<MilanLongInt> &GMate,
                      MilanLongInt *Mate,
                      map<MilanLongInt, MilanLongInt> &Ghost2LocalMap)
{

    /*
#pragma omp critical(Mate)
    {
        if ((node < StartIndex) || (node > EndIndex)) { //Is it a ghost vertex?
            result = GMate[Ghost2LocalMap[node]] >= 0;// Already matched
        } else { //A local vertex
            result = (Mate[node - StartIndex] >= 0); // Already matched
        }

    }
    */
    MilanLongInt val;
    if ((node < StartIndex) || (node > EndIndex)) // if ghost vertex
    {
#pragma omp atomic read
        val = GMate[Ghost2LocalMap[node]];
        return val >= 0; // Already matched
    }

    // If not ghost vertex
#pragma omp atomic read
    val = Mate[node - StartIndex];

    return val >= 0; // Already matched
}
#endif
