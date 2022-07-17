#include "MatchBoxPC.h"

//#define error

void processMatchedVertices(
    MilanLongInt NLVer,
    staticQueue &U,
    staticQueue &privateU,
    MilanLongInt StartIndex,
    MilanLongInt EndIndex,
    MilanLongInt *myCardPtr,
    MilanLongInt *msgIndPtr,
    MilanLongInt *NumMessagesBundledPtr,
    MilanLongInt *SPtr,
    MilanLongInt *verLocPtr,
    MilanLongInt *verLocInd,
    MilanLongInt *verDistance,
    MilanLongInt *PCounter,
    vector<MilanLongInt> &Counter,
    MilanInt myRank,
    MilanInt numProcs,
    MilanLongInt *candidateMate,
    vector<MilanLongInt> &GMate,
    MilanLongInt *Mate,
    map<MilanLongInt, MilanLongInt> &Ghost2LocalMap,
    MilanReal *edgeLocWeight,
    vector<MilanLongInt> &QLocalVtx,
    vector<MilanLongInt> &QGhostVtx,
    vector<MilanLongInt> &QMsgType,
    vector<MilanInt> &QOwner,
    staticQueue &privateQLocalVtx,
    staticQueue &privateQGhostVtx,
    staticQueue &privateQMsgType,
    staticQueue &privateQOwner,
    omp_lock_t *MateLock)
{

    MilanLongInt adj1, adj2, adj11, adj12, k, k1, v = -1, w = -1, ghostOwner;
    int option;
#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << "=========================************===============================" << endl;
    fflush(stdout);
    fflush(stdout);
#endif

#ifdef COUNT_LOCAL_VERTEX
    MilanLongInt localVertices = 0;
#endif
#pragma omp parallel private(k, w, v, k1, adj1, adj2, adj11, adj12, ghostOwner, option) firstprivate(privateU, StartIndex, EndIndex, privateQLocalVtx, privateQGhostVtx, privateQMsgType, privateQOwner) default(shared) num_threads(NUM_THREAD)
    {

        // TODO what would be the optimal UCHUNK
        // TODO refactor
        vector<MilanLongInt> UChunkBeingProcessed;
        UChunkBeingProcessed.reserve(UCHUNK);

        while (!U.empty())
        {

            extractUChunk(UChunkBeingProcessed, U, privateU);

            for (MilanLongInt u : UChunkBeingProcessed)
            {
#ifdef PRINT_DEBUG_INFO_
                cout << "\n(" << myRank << ")u: " << u;
                fflush(stdout);
#endif
                if ((u >= StartIndex) && (u <= EndIndex))
                { // Process Only the Local Vertices

#ifdef COUNT_LOCAL_VERTEX
                    localVertices++;
#endif

                    // Get the Adjacency list for u
                    adj1 = verLocPtr[u - StartIndex]; // Pointer
                    adj2 = verLocPtr[u - StartIndex + 1];
                    for (k = adj1; k < adj2; k++)
                    {
                        option = -1;
                        v = verLocInd[k];

                        if ((v >= StartIndex) && (v <= EndIndex))
                        { // If Local Vertex:

#ifdef PRINT_DEBUG_INFO_
                            cout << "\n(" << myRank << ")v: " << v << " c(v)= " << candidateMate[v - StartIndex] << " Mate[v]: " << Mate[v];
                            fflush(stdout);
#endif

                            // If the current vertex is pointing to a matched vertex and is not matched
                            if (not isAlreadyMatched(v, StartIndex, EndIndex, GMate, Mate, Ghost2LocalMap))
                            {
#pragma omp critical
                                {
                                    if (candidateMate[v - StartIndex] == u)
                                    {
                                        // Start: PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
                                        w = computeCandidateMate(verLocPtr[v - StartIndex],
                                                                 verLocPtr[v - StartIndex + 1],
                                                                 edgeLocWeight, 0,
                                                                 verLocInd,
                                                                 StartIndex,
                                                                 EndIndex,
                                                                 GMate,
                                                                 Mate,
                                                                 Ghost2LocalMap);

                                        candidateMate[v - StartIndex] = w;

#ifdef PRINT_DEBUG_INFO_
                                        cout << "\n(" << myRank << ")" << v << " Points to: " << w;
                                        fflush(stdout);
#endif
                                        // If found a dominating edge:
                                        if (w >= 0)
                                        {

                                            if ((w < StartIndex) || (w > EndIndex))
                                            { // A ghost
#ifdef PRINT_DEBUG_INFO_
                                                cout << "\n(" << myRank << ")Sending a request message:";
                                                cout << "\n(" << myRank << ")Ghost is " << w << " Owner is: " << findOwnerOfGhost(w, verDistance, myRank, numProcs);
#endif
                                                option = 2;

                                                if (candidateMate[NLVer + Ghost2LocalMap[w]] == v)
                                                {
                                                    option = 1;
                                                    Mate[v - StartIndex] = w;     // v is a local vertex
                                                    GMate[Ghost2LocalMap[w]] = v; // w is a ghost vertex

                                                    // Decrement the counter:
                                                    PROCESS_CROSS_EDGE(Counter, Ghost2LocalMap[w], SPtr);
                                                } // End of if CandidateMate[w] = v
                                            }     // End of if a Ghost Vertex
                                            else
                                            { // w is a local vertex
                                                if (candidateMate[w - StartIndex] == v)
                                                {
                                                    option = 3;
                                                    Mate[v - StartIndex] = w; // v is a local vertex
                                                    Mate[w - StartIndex] = v; // w is a local vertex

#ifdef PRINT_DEBUG_INFO_
                                                    cout << "\n(" << myRank << ")MATCH: (" << v << "," << w << ") ";
                                                    fflush(stdout);
#endif
                                                } // End of if(CandidateMate(w) = v
                                            }     // End of Else

                                        } // End of if(w >=0)
                                        else option = 4;// End of Else: w == -1
                                        // End:   PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
                                    }
                                } // End of task
                            }     // End of If (candidateMate[v-StartIndex] == u

                        } // End of if ( (v >= StartIndex) && (v <= EndIndex) ) //If Local Vertex:
                        else
                        { // Neighbor is a ghost vertex

#pragma omp critical
                            {
                                if (candidateMate[NLVer + Ghost2LocalMap[v]] == u)
                                    candidateMate[NLVer + Ghost2LocalMap[v]] = -1;
                                if (v != Mate[u - StartIndex]) option = 5; // u is local
                            } // End of critical
                        }     // End of Else //A Ghost Vertex


                        switch (option)
                        {
                        case -1:
                            // No things to do
                            break;
                        case 1:
                            // Found a dominating edge, it is a ghost and candidateMate[NLVer + Ghost2LocalMap[w]] == v
                            privateU.push_back(v);
                            privateU.push_back(w);
#pragma omp atomic
                            (*myCardPtr)++;
#ifdef PRINT_DEBUG_INFO_
                            cout << "\n(" << myRank << ")MATCH: (" << v << "," << w << ") ";
                            fflush(stdout);
#endif
                        case 2:
                            // Found a dominating edge, it is a ghost
                            ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
                            assert(ghostOwner != -1);
                            assert(ghostOwner != myRank);
#pragma omp atomic
                            PCounter[ghostOwner]++;
#pragma omp atomic
                            (*msgIndPtr)++;
#pragma omp atomic
                            (*NumMessagesBundledPtr)++;
                            privateQLocalVtx.push_back(v);
                            privateQGhostVtx.push_back(w);
                            privateQMsgType.push_back(REQUEST);
                            privateQOwner.push_back(ghostOwner);
                            break;
                        case 3:
                            privateU.push_back(v);
                            privateU.push_back(w);
#pragma omp atomic
                            (*myCardPtr)++;
                            break;
                        case 4:
                            // Could not find a dominating vertex
                            adj11 = verLocPtr[v - StartIndex];
                            adj12 = verLocPtr[v - StartIndex + 1];
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

                                    ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
                                    assert(ghostOwner != -1);
                                    assert(ghostOwner != myRank);
#pragma omp atomic
                                    PCounter[ghostOwner]++;
#pragma omp atomic
                                    (*msgIndPtr)++;
#pragma omp atomic
                                    (*NumMessagesBundledPtr)++;

                                    privateQLocalVtx.push_back(v);
                                    privateQGhostVtx.push_back(w);
                                    privateQMsgType.push_back(FAILURE);
                                    privateQOwner.push_back(ghostOwner);

                                } // End of if(GHOST)
                            }     // End of for loop
                            break;
                        default:
                            
#ifdef PRINT_DEBUG_INFO_
                                    cout << "\n(" << myRank << ")Sending a success message: ";
                                    cout << "\n(" << myRank << ")Ghost is " << v << " Owner is: " << findOwnerOfGhost(v, verDistance, myRank, numProcs) << "\n";
                                    fflush(stdout);
#endif

                                    ghostOwner = findOwnerOfGhost(v, verDistance, myRank, numProcs);
                                    assert(ghostOwner != -1);
                                    assert(ghostOwner != myRank);
#pragma omp atomic
                                    PCounter[ghostOwner]++;
#pragma omp atomic
                                    (*msgIndPtr)++;
#pragma omp atomic
                                    (*NumMessagesBundledPtr)++;
                                    privateQLocalVtx.push_back(u);
                                    privateQGhostVtx.push_back(v);
                                    privateQMsgType.push_back(SUCCESS);
                                    privateQOwner.push_back(ghostOwner);

                            break;
                        } //End of switch

                    } // End of inner for

                    // TODO privateU.size() < UCHUNK could be commented but it generate errors, why?
                    if (privateU.size() > UCHUNK || U.empty())
                    {
#pragma omp critical(U)
                        {
                            while (!privateU.empty())
                                U.push_back(privateU.pop_back());
                        }

#ifndef error
#pragma omp critical(privateMsg)
                        {
                            while (!privateQLocalVtx.empty())
                            {
                                QLocalVtx.push_back(privateQLocalVtx.pop_back());
                                QGhostVtx.push_back(privateQGhostVtx.pop_back());
                                QMsgType.push_back(privateQMsgType.pop_back());
                                QOwner.push_back(privateQOwner.pop_back());
                            }
                        }

#endif
                    } // End of private.size()
                }
            } // End of outer for
        }     // End of while ( !U.empty() )
        queuesTransfer(U, privateU, QLocalVtx,
                       QGhostVtx,
                       QMsgType, QOwner, privateQLocalVtx,
                       privateQGhostVtx,
                       privateQMsgType,
                       privateQOwner);

#ifdef COUNT_LOCAL_VERTEX
        printf("Count local vertexes: %ld for thread %d of processor %d\n",
               localVertices,
               omp_get_thread_num(),
               myRank);

#endif
    } // End of parallel region
}
