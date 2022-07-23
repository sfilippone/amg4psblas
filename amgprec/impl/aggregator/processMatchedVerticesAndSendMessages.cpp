#include "MatchBoxPC.h"

void processMatchedVerticesAndSendMessages(
    MilanLongInt NLVer,
    vector<MilanLongInt> &UChunkBeingProcessed,
    staticQueue &U,
    staticQueue &privateU,
    MilanLongInt StartIndex,
    MilanLongInt EndIndex,
    MilanLongInt *myCard,
    MilanLongInt *msgInd,
    MilanLongInt *NumMessagesBundled,
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
    vector<MilanLongInt> &privateQLocalVtx,
    vector<MilanLongInt> &privateQGhostVtx,
    vector<MilanLongInt> &privateQMsgType,
    vector<MilanInt> &privateQOwner,
    MPI_Comm comm,
    MilanLongInt *msgActual,
    vector<MilanLongInt> &Message)
{

    MilanLongInt initialSize = QLocalVtx.size();
    MilanLongInt adj1, adj2, adj11, adj12, k, k1, v = -1, w = -1, ghostOwner;
    int option;
    MilanLongInt mateVal;

#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << "=========================************===============================" << endl;
    fflush(stdout);
    fflush(stdout);
#endif

#ifdef COUNT_LOCAL_VERTEX
    MilanLongInt localVertices = 0;
#endif
#pragma omp parallel private(k, w, v, k1, adj1, adj2, adj11, adj12, ghostOwner, option)                                                                             \
    firstprivate(Message, privateU, StartIndex, EndIndex, privateQLocalVtx, privateQGhostVtx, privateQMsgType, privateQOwner, UChunkBeingProcessed) default(shared) \
        num_threads(NUM_THREAD)                                                                                                                                     \
            reduction(+                                                                                                                                             \
                      : msgInd[:1], PCounter                                                                                                                        \
                      [:numProcs], myCard                                                                                                                           \
                      [:1], NumMessagesBundled                                                                                                                      \
                      [:1], msgActual                                                                                                                               \
                      [:1])
    {

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
#pragma omp atomic read
                            mateVal = Mate[v - StartIndex];
                            // If the current vertex is pointing to a matched vertex and is not matched
                            if (mateVal < 0)
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
                                        }         // End of if(w >=0)
                                        else
                                            option = 4; // End of Else: w == -1
                                        // End:   PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
                                    } // End of If (candidateMate[v-StartIndex] == u
                                }     // End of task
                            }         // mateval < 0
                        }             // End of if ( (v >= StartIndex) && (v <= EndIndex) ) //If Local Vertex:
                        else
                        { // Neighbor is a ghost vertex

#pragma omp critical
                            {
                                if (candidateMate[NLVer + Ghost2LocalMap[v]] == u)
                                    candidateMate[NLVer + Ghost2LocalMap[v]] = -1;
                                if (v != Mate[u - StartIndex])
                                    option = 5; // u is local
                            }                   // End of critical
                        }                       // End of Else //A Ghost Vertex

                        switch (option)
                        {
                        case -1:
                            // No things to do
                            break;
                        case 1:
                            // Found a dominating edge, it is a ghost and candidateMate[NLVer + Ghost2LocalMap[w]] == v
                            privateU.push_back(v);
                            privateU.push_back(w);
                            (*myCard)++;
#ifdef PRINT_DEBUG_INFO_
                            cout << "\n(" << myRank << ")MATCH: (" << v << "," << w << ") ";
                            fflush(stdout);
#endif
                            // Decrement the counter:
                            PROCESS_CROSS_EDGE(&Counter[Ghost2LocalMap[w]], SPtr);
                        case 2:

                            // Found a dominating edge, it is a ghost
                            ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);

                            // Build the Message Packet:
                            // Message[0] = v;       // LOCAL
                            // Message[1] = w;       // GHOST
                            // Message[2] = REQUEST; // TYPE
                            // Send a Request (Asynchronous)
                            // MPI_Bsend(&Message[0], 3, TypeMap<MilanLongInt>(), ghostOwner, ComputeTag, comm);

                            (*msgActual)++;
                            (*msgInd)++;

                            privateQLocalVtx.push_back(v);
                            privateQGhostVtx.push_back(w);
                            privateQMsgType.push_back(REQUEST);
                            privateQOwner.push_back(ghostOwner);
                            break;
                        case 3:
                            privateU.push_back(v);
                            privateU.push_back(w);
                            (*myCard)++;
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

                                    // Build the Message Packet:
                                    // Message[0] = v;       // LOCAL
                                    // Message[1] = w;       // GHOST
                                    // Message[2] = FAILURE; // TYPE
                                    // Send a Request (Asynchronous)
                                    // MPI_Bsend(&Message[0], 3, TypeMap<MilanLongInt>(), ghostOwner, ComputeTag, comm);

                                    (*msgActual)++;
                                    (*msgInd)++;

                                    privateQLocalVtx.push_back(v);
                                    privateQGhostVtx.push_back(w);
                                    privateQMsgType.push_back(FAILURE);
                                    privateQOwner.push_back(ghostOwner);

                                } // End of if(GHOST)
                            }     // End of for loop
                            break;
                        case 5:
                        default:

#ifdef PRINT_DEBUG_INFO_
                            cout << "\n(" << myRank << ")Sending a success message: ";
                            cout << "\n(" << myRank << ")Ghost is " << v << " Owner is: " << findOwnerOfGhost(v, verDistance, myRank, numProcs) << "\n";
                            fflush(stdout);
#endif

                            ghostOwner = findOwnerOfGhost(v, verDistance, myRank, numProcs);

                            // Build the Message Packet:
                            // Message[0] = u;       // LOCAL
                            // Message[1] = v;       // GHOST
                            // Message[2] = SUCCESS; // TYPE
                            // Send a Request (Asynchronous)
                            // MPI_Bsend(&Message[0], 3, TypeMap<MilanLongInt>(), ghostOwner, ComputeTag, comm);

                            (*msgActual)++;
                            (*msgInd)++;

                            privateQLocalVtx.push_back(u);
                            privateQGhostVtx.push_back(v);
                            privateQMsgType.push_back(SUCCESS);
                            privateQOwner.push_back(ghostOwner);

                            break;
                        } // End of switch
                    }     // End of inner for
                }
            } // End of outer for

            queuesTransfer(U, privateU, QLocalVtx,
                           QGhostVtx,
                           QMsgType, QOwner, privateQLocalVtx,
                           privateQGhostVtx,
                           privateQMsgType,
                           privateQOwner);

        } // End of while ( !U.empty() )

#ifdef COUNT_LOCAL_VERTEX
        printf("Count local vertexes: %ld for thread %d of processor %d\n",
               localVertices,
               omp_get_thread_num(),
               myRank);

#endif
    } // End of parallel region

    // Send the messages
    for (int i = initialSize; i < QOwner.size(); i++)
    {

        Message[0] = QLocalVtx[i];
        Message[1] = QGhostVtx[i];
        Message[2] = QMsgType[i];
        ghostOwner = QOwner[i];

        MPI_Bsend(&Message[0], 3, TypeMap<MilanLongInt>(), ghostOwner, ComputeTag, comm);
    }
}
