#include "MatchBoxPC.h"

/**
 * Execute the research fr the Candidate Mate without controlling if the vertices are already matched.
 * Returns the vertices with the highest weight
 * @param adj1
 * @param adj2
 * @param verLocInd
 * @param edgeLocWeight
 * @return
 */
inline MilanLongInt firstComputeCandidateMate(MilanLongInt adj1,
                                              MilanLongInt adj2,
                                              MilanLongInt* verLocInd,
                                              MilanReal* edgeLocWeight)
{
    MilanInt w = -1;
    MilanReal heaviestEdgeWt = MilanRealMin; //Assign the smallest Value possible first LDBL_MIN
    int finalK;
    for (int k = adj1; k < adj2; k++) {

        if ((edgeLocWeight[k] > heaviestEdgeWt) ||
            ((edgeLocWeight[k] == heaviestEdgeWt) && (w < verLocInd[k]))) {
            heaviestEdgeWt = edgeLocWeight[k];
            w = verLocInd[k];
            finalK = k;
        }
    } //End of for loop
    return finalK;
}

/**
 * //TODO documentation
 * @param adj1
 * @param adj2
 * @param edgeLocWeight
 * @param k
 * @param verLocInd
 * @param StartIndex
 * @param EndIndex
 * @param GMate
 * @param Mate
 * @param Ghost2LocalMap
 * @return
 */
inline MilanLongInt computeCandidateMate(MilanLongInt adj1,
                                         MilanLongInt adj2,
                                         MilanReal* edgeLocWeight,
                                         MilanLongInt k,
                                         MilanLongInt* verLocInd,
                                         MilanLongInt StartIndex,
                                         MilanLongInt EndIndex,
                                         vector <MilanLongInt>& GMate,
                                         MilanLongInt* Mate,
                                         map <MilanLongInt, MilanLongInt>& Ghost2LocalMap)
{
    MilanInt w = -1;
    MilanReal heaviestEdgeWt = MilanRealMin; //Assign the smallest Value possible first LDBL_MIN
    for (k = adj1; k < adj2; k++) {
        if (isAlreadyMatched(verLocInd[k], StartIndex, EndIndex, GMate, Mate, Ghost2LocalMap)) continue;

        if ((edgeLocWeight[k] > heaviestEdgeWt) ||
            ((edgeLocWeight[k] == heaviestEdgeWt) && (w < verLocInd[k]))) {
            heaviestEdgeWt = edgeLocWeight[k];
            w = verLocInd[k];
        }
    } //End of for loop
    return w;
}