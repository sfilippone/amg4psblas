#include "MatchBoxPC.h"

void processMessages(
    MilanLongInt NLVer,
    MilanLongInt *Mate,
    MilanLongInt *candidateMate,
    map<MilanLongInt, MilanLongInt> &Ghost2LocalMap,
    vector<MilanLongInt> &GMate,
    vector<MilanLongInt> &Counter,
    MilanLongInt StartIndex,
    MilanLongInt EndIndex,
    MilanLongInt *myCardPtr,
    MilanLongInt *msgIndPtr,
    MilanLongInt *msgActualPtr,
    MilanReal *edgeLocWeight,
    MilanLongInt *verDistance,
    MilanLongInt *verLocPtr,
    MilanLongInt k,
    MilanLongInt *verLocInd,
    MilanInt numProcs,
    MilanInt myRank,
    MPI_Comm comm,
    vector<MilanLongInt> &Message,
    MilanLongInt numGhostEdges,
    MilanLongInt u,
    MilanLongInt v,
    MilanLongInt *S,
    staticQueue &U)
{

    MilanInt Sender;
    MPI_Status computeStatus;
    MilanLongInt bundleSize, myCard = *myCardPtr, msgInd = *msgIndPtr, msgActual = *msgActualPtr, w;
    MilanLongInt adj11, adj12, k1;
    MilanLongInt ghostOwner;
    int error_codeC;
    error_codeC = MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    char error_message[MPI_MAX_ERROR_STRING];
    int message_length;
    MilanLongInt message_type = 0;

    // Buffer to receive bundled messages
    // Maximum messages that can be received from any processor is
    // twice the edge cut: REQUEST; REQUEST+(FAILURE/SUCCESS)
    vector<MilanLongInt> ReceiveBuffer;
    try
    {
        ReceiveBuffer.reserve(numGhostEdges * 2 * 3); // Three integers per cross edge
    }
    catch (length_error)
    {
        cout << "Error in function algoDistEdgeApproxDominatingEdgesMessageBundling: \n";
        cout << "Not enough memory to allocate the internal variables \n";
        exit(1);
    }

#ifdef PRINT_DEBUG_INFO_
    cout
        << "\n(" << myRank << "=========================************===============================" << endl;
    fflush(stdout);
    fflush(stdout);
#endif
#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << ")About to begin Message processing phase ... *S=" << *S << endl;
    fflush(stdout);
#endif
#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << "=========================************===============================" << endl;
    fflush(stdout);
    fflush(stdout);
#endif
    // BLOCKING RECEIVE:
#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << " Waiting for blocking receive..." << endl;
    fflush(stdout);
    fflush(stdout);
#endif

    error_codeC = MPI_Recv(&Message[0], 3, TypeMap<MilanLongInt>(), MPI_ANY_SOURCE, ComputeTag, comm, &computeStatus);
    if (error_codeC != MPI_SUCCESS)
    {
        MPI_Error_string(error_codeC, error_message, &message_length);
        cout << "\n*Error in call to MPI_Receive on Slave: " << error_message << "\n";
        fflush(stdout);
    }
    Sender = computeStatus.MPI_SOURCE;

#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << ")Received message from Process " << Sender << " Type= " << Message[2] << endl;
    fflush(stdout);
#endif

    if (Message[2] == SIZEINFO)
    {

#ifdef PRINT_DEBUG_INFO_
        cout << "\n(" << myRank << ")Received bundled message from Process " << Sender << " Size= " << Message[0] << endl;
        fflush(stdout);
#endif
        bundleSize = Message[0]; //#of integers in the message
        // Build the Message Buffer:
        if (!ReceiveBuffer.empty())
            ReceiveBuffer.clear();            // Empty it out first
        ReceiveBuffer.resize(bundleSize, -1); // Initialize
#ifdef PRINT_DEBUG_INFO_
        cout << "\n(" << myRank << ")Message Bundle Before: " << endl;
        for (i = 0; i < bundleSize; i++)
            cout << ReceiveBuffer[i] << ",";
        cout << endl;
        fflush(stdout);
#endif
        // Receive the message
        error_codeC = MPI_Recv(&ReceiveBuffer[0], bundleSize, TypeMap<MilanLongInt>(), Sender, BundleTag, comm, &computeStatus);
        if (error_codeC != MPI_SUCCESS)
        {
            MPI_Error_string(error_codeC, error_message, &message_length);
            cout << "\n*Error in call to MPI_Receive on processor " << myRank << " Error: " << error_message << "\n";
            fflush(stdout);
        }
#ifdef PRINT_DEBUG_INFO_
        cout << "\n(" << myRank << ")Message Bundle After: " << endl;
        for (i = 0; i < bundleSize; i++)
            cout << ReceiveBuffer[i] << ",";
        cout << endl;
        fflush(stdout);
#endif
    }
    else
    { // Just a single message:
#ifdef PRINT_DEBUG_INFO_
        cout << "\n(" << myRank << ")Received regular message from Process " << Sender << " u= " << Message[0] << " v= " << Message[1] << endl;
        fflush(stdout);
#endif
        // Add the current message to Queue:
        bundleSize = 3; //#of integers in the message
        // Build the Message Buffer:
        if (!ReceiveBuffer.empty())
            ReceiveBuffer.clear();            // Empty it out first
        ReceiveBuffer.resize(bundleSize, -1); // Initialize

        ReceiveBuffer[0] = Message[0]; // u
        ReceiveBuffer[1] = Message[1]; // v
        ReceiveBuffer[2] = Message[2]; // message_type
    }

#ifdef DEBUG_GHOST_
    if ((v < StartIndex) || (v > EndIndex))
    {
        cout << "\n(" << myRank << ") From ReceiveBuffer: This should not happen: u= " << u << " v= " << v << " Type= " << message_type << " StartIndex " << StartIndex << " EndIndex " << EndIndex << endl;
        fflush(stdout);
    }
#endif
#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << ")Processing message: u= " << u << " v= " << v << " Type= " << message_type << endl;
    fflush(stdout);
#endif


    //Most of the time bundleSize == 3, thus, it's not worth parallelizing thi loop
    for (MilanLongInt bundleCounter = 3; bundleCounter < bundleSize + 3; bundleCounter += 3)
    {
        u = ReceiveBuffer[bundleCounter - 3];            // GHOST
        v = ReceiveBuffer[bundleCounter - 2];            // LOCAL
        message_type = ReceiveBuffer[bundleCounter - 1]; // TYPE

        // CASE I: REQUEST
        if (message_type == REQUEST)
        {
#ifdef PRINT_DEBUG_INFO_
            cout << "\n(" << myRank << ")Message type is REQUEST" << endl;
            fflush(stdout);
#endif
#ifdef DEBUG_GHOST_
            if ((v < 0) || (v < StartIndex) || ((v - StartIndex) > NLVer))
            {
                cout << "\n(" << myRank << ") case 1 Bad address " << v << " " << StartIndex << " " << v - StartIndex << " " << NLVer << endl;
                fflush(stdout);
            }

#endif

            if (Mate[v - StartIndex] == -1)
            {                                                 // Process only if not already matched  (v is local)
                candidateMate[NLVer + Ghost2LocalMap[u]] = v; // Set CandidateMate for the ghost
                if (candidateMate[v - StartIndex] == u)
                {
                    GMate[Ghost2LocalMap[u]] = v; // u is ghost
                    Mate[v - StartIndex] = u;     // v is local
                    U.push_back(v);
                    U.push_back(u);
                    myCard++;
#ifdef PRINT_DEBUG_INFO_
                    cout << "\n(" << myRank << ")MATCH: (" << v << "," << u << ") " << endl;
                    fflush(stdout);
#endif

                    PROCESS_CROSS_EDGE(&Counter[Ghost2LocalMap[u]], S);
                } // End of if ( candidateMate[v-StartIndex] == u )e
            }     // End of if ( Mate[v] == -1 )
        }         // End of REQUEST
        else
        { // CASE II: SUCCESS
            if (message_type == SUCCESS)
            {
#ifdef PRINT_DEBUG_INFO_
                cout << "\n(" << myRank << ")Message type is SUCCESS" << endl;
                fflush(stdout);
#endif
                GMate[Ghost2LocalMap[u]] = EndIndex + 1; // Set a Dummy Mate to make sure that we do not (u is a ghost) process it again
                PROCESS_CROSS_EDGE(&Counter[Ghost2LocalMap[u]], S);
#ifdef DEBUG_GHOST_
                if ((v < 0) || (v < StartIndex) || ((v - StartIndex) > NLVer))
                {
                    cout << "\n(" << myRank << ") case 2  Bad address " << v << " " << StartIndex << " " << v - StartIndex << " " << NLVer << endl;
                    fflush(stdout);
                }
#endif
                if (Mate[v - StartIndex] == -1)
                { // Process only if not already matched ( v is local)
                    if (candidateMate[v - StartIndex] == u)
                    {
                        // Start: PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
                        w = computeCandidateMate(verLocPtr[v - StartIndex], verLocPtr[v - StartIndex + 1], edgeLocWeight, k, verLocInd, StartIndex, EndIndex, GMate, Mate, Ghost2LocalMap);
                        candidateMate[v - StartIndex] = w;
#ifdef PRINT_DEBUG_INFO_
                        cout << "\n(" << myRank << ")" << v << " Points to: " << w << endl;
                        fflush(stdout);
#endif
                        // If found a dominating edge:
                        if (w >= 0)
                        {
                            if ((w < StartIndex) || (w > EndIndex))
                            { // w is a ghost
                                // Build the Message Packet:
                                Message[0] = v;       // LOCAL
                                Message[1] = w;       // GHOST
                                Message[2] = REQUEST; // TYPE
                                                      // Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
                                cout << "\n(" << myRank << ")Sending a request message: ";
                                cout << "\n(" << myRank << ")Ghost is " << w << " Owner is: " << findOwnerOfGhost(w, verDistance, myRank, numProcs) << endl;
                                fflush(stdout);
#endif
                                ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
                                assert(ghostOwner != -1);
                                assert(ghostOwner != myRank);

                                MPI_Bsend(&Message[0], 3, TypeMap<MilanLongInt>(), ghostOwner, ComputeTag, comm);
                                msgInd++;
                                msgActual++;
                                if (candidateMate[NLVer + Ghost2LocalMap[w]] == v)
                                {
                                    Mate[v - StartIndex] = w;     // v is local
                                    GMate[Ghost2LocalMap[w]] = v; // w is ghost
                                    U.push_back(v);
                                    U.push_back(w);
                                    myCard++;
#ifdef PRINT_DEBUG_INFO_
                                    cout << "\n(" << myRank << ")MATCH: (" << v << "," << w << ") " << endl;
                                    fflush(stdout);
#endif

                                    PROCESS_CROSS_EDGE(&Counter[Ghost2LocalMap[w]], S);
                                } // End of if CandidateMate[w] = v
                            }     // End of if a Ghost Vertex
                            else
                            { // w is a local vertex
                                if (candidateMate[w - StartIndex] == v)
                                {
                                    Mate[v - StartIndex] = w; // v is local
                                    Mate[w - StartIndex] = v; // w is local
                                    // Q.push_back(u);
                                    U.push_back(v);
                                    U.push_back(w);
                                    myCard++;
#ifdef PRINT_DEBUG_INFO_
                                    cout << "\n(" << myRank << ")MATCH: (" << v << "," << w << ") " << endl;
                                    fflush(stdout);
#endif
                                } // End of if(CandidateMate(w) = v
                            }     // End of Else
                        }         // End of if(w >=0)
                        else
                        { // No dominant edge found
                            adj11 = verLocPtr[v - StartIndex];
                            adj12 = verLocPtr[v - StartIndex + 1];
                            for (k1 = adj11; k1 < adj12; k1++)
                            {
                                w = verLocInd[k1];
                                if ((w < StartIndex) || (w > EndIndex))
                                { // A ghost
                                    // Build the Message Packet:
                                    Message[0] = v;       // LOCAL
                                    Message[1] = w;       // GHOST
                                    Message[2] = FAILURE; // TYPE
                                                          // Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
                                    cout << "\n(" << myRank << ")Sending a failure message: ";
                                    cout << "\n(" << myRank << ")Ghost is " << w << " Owner is: " << findOwnerOfGhost(w, verDistance, myRank, numProcs) << endl;
                                    fflush(stdout);
#endif
                                    ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
                                    assert(ghostOwner != -1);
                                    assert(ghostOwner != myRank);
                                    MPI_Bsend(&Message[0], 3, TypeMap<MilanLongInt>(), ghostOwner, ComputeTag, comm);
                                    msgInd++;
                                    msgActual++;
                                } // End of if(GHOST)
                            }     // End of for loop
                        }         // End of Else: w == -1
                        // End:   PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
                    } // End of if ( candidateMate[v-StartIndex] == u )
                }     // End of if ( Mate[v] == -1 )
            }         // End of if ( message_type == SUCCESS )
            else
            { // CASE III: FAILURE
#ifdef PRINT_DEBUG_INFO_
                cout << "\n(" << myRank << ")Message type is FAILURE" << endl;
                fflush(stdout);
#endif
                GMate[Ghost2LocalMap[u]] = EndIndex + 1; // Set a Dummy Mate to make sure that we do not (u is a ghost) process this anymore
                PROCESS_CROSS_EDGE(&Counter[Ghost2LocalMap[u]], S); // Decrease the counter
            }                                                      // End of else: CASE III
        }                                                          // End of else: CASE I
    }

    *myCardPtr = myCard;
    *msgIndPtr = msgInd;
    *msgActualPtr = msgActual;
    return;
}