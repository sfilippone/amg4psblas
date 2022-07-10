#include "MatchBoxPC.h"
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <map>
#include <vector>
#include "primitiveDataTypeDefinitions.h"
#include "dataStrStaticQueue.h"
#include "omp.h"
#include "queueTransfer.cpp"
#include "processCrossEdge.cpp"

inline void PARALLEL_PROCESS_EXPOSED_VERTEX_B(MilanLongInt NLVer,
                                              MilanLongInt *candidateMate,
                                              MilanLongInt *verLocInd,
                                              MilanLongInt *verLocPtr,
                                              MilanLongInt StartIndex,
                                              MilanLongInt EndIndex,
                                              MilanLongInt *Mate,
                                              vector<MilanLongInt> &GMate,
                                              map<MilanLongInt, MilanLongInt> &Ghost2LocalMap,
                                              MilanReal *edgeLocWeight,
                                              MilanLongInt *myCardPtr,
                                              MilanLongInt *msgIndPtr,
                                              MilanLongInt *NumMessagesBundledPtr,
                                              MilanLongInt *SPtr,
                                              MilanLongInt *verDistance,
                                              MilanLongInt *PCounter,
                                              vector<MilanLongInt> &Counter,
                                              MilanInt myRank,
                                              MilanInt numProcs,
                                              staticQueue &U,
                                              staticQueue &privateU,
                                              vector<MilanLongInt> &QLocalVtx,
                                              vector<MilanLongInt> &QGhostVtx,
                                              vector<MilanLongInt> &QMsgType,
                                              vector<MilanInt> &QOwner,
                                              staticQueue &privateQLocalVtx,
                                              staticQueue &privateQGhostVtx,
                                              staticQueue &privateQMsgType,
                                              staticQueue &privateQOwner)
{

    //TODO define all the constants in a single place!
    const MilanLongInt REQUEST = 1;
    const MilanLongInt SUCCESS = 2;
    const MilanLongInt FAILURE = 3;
    const MilanLongInt SIZEINFO = 4;
    MilanLongInt v = -1, k = -1, w = -1, adj11 = 0, adj12 = 0, k1 = 0, S = *SPtr;
    MilanLongInt myCard = 0, msgInd = 0;
    MilanLongInt NumMessagesBundled = 0;
    MilanInt ghostOwner = 0;

#pragma omp parallel private(k, w, v, k1, adj11, adj12, ghostOwner) firstprivate(privateU, StartIndex, EndIndex, privateQLocalVtx, privateQGhostVtx, privateQMsgType, privateQOwner) default(shared) num_threads(4)
    {
#pragma omp for reduction(+ \
                          : msgInd, NumMessagesBundled, myCard, PCounter[:numProcs]) schedule(static)
        for (v = 0; v < NLVer; v++)
        {
            // Start: PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
            k = candidateMate[v];
            candidateMate[v] = verLocInd[k];
            w = candidateMate[v];

#ifdef PRINT_DEBUG_INFO_
            cout << "\n(" << myRank << ")Processing: " << v + StartIndex << endl;
            fflush(stdout);
#endif

#ifdef PRINT_DEBUG_INFO_
            cout << "\n(" << myRank << ")" << v + StartIndex << " Points to: " << w;
            fflush(stdout);
#endif
            // If found a dominating edge:
            if (w >= 0)
            {

                if (isAlreadyMatched(verLocInd[k], StartIndex, EndIndex, GMate, Mate, Ghost2LocalMap))
                {
                    w = computeCandidateMate(verLocPtr[v],
                                             verLocPtr[v + 1],
                                             edgeLocWeight, 0,
                                             verLocInd,
                                             StartIndex,
                                             EndIndex,
                                             GMate,
                                             Mate,
                                             Ghost2LocalMap);
                    candidateMate[v] = w;
                }

                if (w >= 0)
                {

                    myCard++;
                    if ((w < StartIndex) || (w > EndIndex))
                    { // w is a ghost vertex
#ifdef PRINT_DEBUG_INFO_
                        cout << "\n(" << myRank << ")Sending a request message (291):";
                        cout << "\n(" << myRank << ")Local is: " << v + StartIndex << " Ghost is " << w << " Owner is: " << findOwnerOfGhost(w, verDistance, myRank, numProcs) << endl;
                        fflush(stdout);
#endif

                        msgInd++;
                        NumMessagesBundled++;
                        ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
                        assert(ghostOwner != -1);
                        assert(ghostOwner != myRank);
                        PCounter[ghostOwner]++;

                        privateQLocalVtx.push_back(v + StartIndex);
                        privateQGhostVtx.push_back(w);
                        privateQMsgType.push_back(REQUEST);
                        privateQOwner.push_back(ghostOwner);
                        

                        if (candidateMate[NLVer + Ghost2LocalMap[w]] == v + StartIndex)
                        {

                            privateU.push_back(v + StartIndex);
                            privateU.push_back(w);
                            Mate[v] = w;
                            // FIXME could this instruction create errors?
                            GMate[Ghost2LocalMap[w]] = v + StartIndex; // w is a Ghost

#ifdef PRINT_DEBUG_INFO_
                            cout << "\n(" << myRank << ")MATCH: (" << v + StartIndex << "," << w << ")";
                            fflush(stdout);
#endif

                            //TODO refactor this!!
                            // Decrement the counter:
                            PROCESS_CROSS_EDGE(Counter, Ghost2LocalMap, w, &S);
                        } // End of if CandidateMate[w] = v

                    } // End of if a Ghost Vertex
                    else
                    { // w is a local vertex

                        if (candidateMate[w - StartIndex] == (v + StartIndex))
                        {
                            privateU.push_back(v + StartIndex);
                            privateU.push_back(w);

                            Mate[v] = w; // v is local
                            // FIXME this instruction could create errors
                            Mate[w - StartIndex] = v + StartIndex; // w is local

#ifdef PRINT_DEBUG_INFO_
                            cout << "\n(" << myRank << ")MATCH: (" << v + StartIndex << "," << w << ") ";
                            fflush(stdout);
#endif

                        } // End of if ( candidateMate[w-StartIndex] == (v+StartIndex) )
                    }     // End of Else

                    continue;
                } // End of second if

            } // End of if(w >=0)

            // This piece of code is executed a really small amount of times, I will not allocate a
            // huge amount of memory for the private data structures.
            adj11 = verLocPtr[v];
            adj12 = verLocPtr[v + 1];
            for (k1 = adj11; k1 < adj12; k1++)
            {
                w = verLocInd[k1];
                if ((w < StartIndex) || (w > EndIndex))
                { // A ghost

#ifdef PRINT_DEBUG_INFO_
                    cout << "\n(" << myRank << ")Sending a failure message: ";
                    cout << "\n(" << myRank << ")Ghost is " << w << " Owner is: " << findOwnerOfGhost(w, verDistance, myRank, numProcs);
                    fflush(stdout);
#endif

                    msgInd++;
                    NumMessagesBundled++;
                    ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
                    assert(ghostOwner != -1);
                    assert(ghostOwner != myRank);
                    PCounter[ghostOwner]++;

                    privateQLocalVtx.push_back(v + StartIndex);
                    privateQGhostVtx.push_back(w);
                    privateQMsgType.push_back(FAILURE);
                    privateQOwner.push_back(ghostOwner);

                } // End of if(GHOST)
            }     // End of for loop
            // End:   PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
        } // End of for ( v=0; v < NLVer; v++ )

        queuesTransfer(U, privateU, QLocalVtx,
                       QGhostVtx,
                       QMsgType, QOwner, privateQLocalVtx,
                       privateQGhostVtx,
                       privateQMsgType,
                       privateQOwner);

//TODO move this outside of the parallel region!!
#pragma omp master
        {
            *myCardPtr = myCard;
            *msgIndPtr = msgInd;
            *NumMessagesBundledPtr = NumMessagesBundled;
            *SPtr = S;
        }

    } // End of parallel region
}