#include "MatchBoxPC.h"

// ***********************************************************************
//
//        MatchboxP: A C++ library for approximate weighted matching
//               Mahantesh Halappanavar (hala@pnnl.gov)
//               Pacific Northwest National Laboratory
//
// ***********************************************************************
//
//       Copyright (2021) Battelle Memorial Institute
//                      All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
//////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// DOMINATING EDGES MODEL ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
/* Function	: algoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMate()
 *
 * Date     : New update: Feb 17, 2019, Richland, Washington.
 * Date		: Original development: May 17, 2009, E&CS Bldg.
 *
 * Purpose	: Compute Approximate Maximum Weight Matching in Linear Time
 *
 * Args		: inputMatrix - instance of Compressed-Col format of Matrix
 *                Mate - The Mate array
 *
 * Returns	: By Value: (void)
 *            By Reference: Mate
 *
 * Comments	: 1/2 Approx Algorithm. Picks the locally available heaviest edge.
 *                Assumption: The Mate Array is empty.
 */

/*
 NLVer = #of vertices, NLEdge = #of edges
 CSR/CSC/Compressed format: verLocPtr = Pointer, verLocInd = Index, edgeLocWeight = edge weights (positive real numbers)
 verDistance = A vector of size |P|+1 containing the cumulative number of vertices per process
 Mate = A vector of size |V_p| (local subgraph) to store the output (matching)
 MPI: myRank, numProcs, comm,
 Statistics: msgIndSent, msgActualSent, msgPercent : Size: |P| number of processes in the comm-world
 Statistics: ph0_time, ph1_time, ph2_time: Runtimes
 Statistics: ph1_card, ph2_card : Size: |P| number of processes in the comm-world (number of matched edges in Phase 1 and Phase 2)
 */

#ifdef SERIAL_MPI
#else

// DOUBLE PRECISION VERSION
// WARNING: The vertex block on a given rank is contiguous
void dalgoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMateCMP(
    MilanLongInt NLVer, MilanLongInt NLEdge,
    MilanLongInt *verLocPtr, MilanLongInt *verLocInd,
    MilanReal *edgeLocWeight,
    MilanLongInt *verDistance,
    MilanLongInt *Mate,
    MilanInt myRank, MilanInt numProcs, MPI_Comm comm,
    MilanLongInt *msgIndSent, MilanLongInt *msgActualSent,
    MilanReal *msgPercent,
    MilanReal *ph0_time, MilanReal *ph1_time, MilanReal *ph2_time,
    MilanLongInt *ph1_card, MilanLongInt *ph2_card)
{

    /*
     * verDistance: it's a vector long as the number of processors.
     *              verDistance[i] contains the first node index of the i-th processor
     *              verDistance[i + 1] contains the last node index of the i-th processor
     * NLVer: number of elements in the LocPtr
     * NLEdge: number of edges assigned to the current processor
     *
     * Contains the portion of matrix assigned to the processor in
     * Yale notation
     * verLocInd: contains the positions on row of the matrix
     * verLocPtr: i-th value is the position of the first element on the i-th row and
     *            i+1-th value is the position of the first element on the i+1-th row
     */

#if !defined(SERIAL_MPI)
#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << ")Within algoEdgeApproxDominatingEdgesLinearSearchMessageBundling()";
    fflush(stdout);
#endif

#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << ") verDistance [" << verDistance[0] << "," << verDistance[1] << "," << verDistance[2] << "," << verDistance[3] << "]";
    fflush(stdout);
#endif
#ifdef DEBUG_HANG_
    if (myRank == 0)
        cout << "\n(" << myRank << ") verDistance [" << verDistance[0] << "," << verDistance[1] << "," << verDistance[2] << "," << verDistance[3] << "]";
    fflush(stdout);
#endif

    MilanLongInt StartIndex = verDistance[myRank];       // The starting vertex owned by the current rank
    MilanLongInt EndIndex = verDistance[myRank + 1] - 1; // The ending vertex owned by the current rank

    MPI_Status computeStatus;

    MilanLongInt msgActual = 0, msgInd = 0;
    MilanReal heaviestEdgeWt = 0.0f; // Assumes positive weight
    MilanReal startTime, finishTime;

    startTime = MPI_Wtime();

    // Data structures for sending and receiving messages:
    vector<MilanLongInt> Message; // [ u, v, message_type ]
    Message.resize(3, -1);
    // Data structures for Message Bundling:
    // Although up to two messages can be sent along any cross edge,
    // only one message will be sent in the initialization phase -
    // one of: REQUEST/FAILURE/SUCCESS
    vector<MilanLongInt> QLocalVtx, QGhostVtx, QMsgType;
    vector<MilanInt> QOwner; // Changed by Fabio to be an integer, addresses needs to be integers!

    MilanLongInt *PCounter = new MilanLongInt[numProcs];
    for (int i = 0; i < numProcs; i++)
        PCounter[i] = 0;

    MilanLongInt NumMessagesBundled = 0;
    // TODO when the last computational section will be refactored this could be eliminated
    MilanInt ghostOwner = 0; // Changed by Fabio to be an integer, addresses needs to be integers!
    MilanLongInt *candidateMate = nullptr;
#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << ")NV: " << NLVer << "  Edges: " << NLEdge;
    fflush(stdout);
    cout << "\n(" << myRank << ")StartIndex: " << StartIndex << "  EndIndex: " << EndIndex;
    fflush(stdout);
#endif
    // Other Variables:
    MilanLongInt u = -1, v = -1, w = -1, i = 0;
    MilanLongInt k = -1, adj1 = -1, adj2 = -1;
    MilanLongInt k1 = -1, adj11 = -1, adj12 = -1;
    MilanLongInt myCard = 0;

    // Build the Ghost Vertex Set: Vg
    map<MilanLongInt, MilanLongInt> Ghost2LocalMap;       // Map each ghost vertex to a local vertex
    vector<MilanLongInt> Counter;                         // Store the edge count for each ghost vertex
    MilanLongInt numGhostVertices = 0, numGhostEdges = 0; // Number of Ghost vertices

#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << ")About to compute Ghost Vertices...";
    fflush(stdout);
#endif
#ifdef DEBUG_HANG_
    if (myRank == 0)
        cout << "\n(" << myRank << ")About to compute Ghost Vertices...";
    fflush(stdout);
#endif

    // Define Adjacency Lists for Ghost Vertices:
    // cout<<"Building Ghost data structures ... \n\n";
    vector<MilanLongInt> verGhostPtr, verGhostInd, tempCounter;
    // Mate array for ghost vertices:
    vector<MilanLongInt> GMate; // Proportional to the number of ghost vertices
    MilanLongInt S;
    MilanLongInt privateMyCard = 0;
    staticQueue U, privateU, privateQLocalVtx, privateQGhostVtx, privateQMsgType, privateQOwner;
    vector<MilanLongInt> PCumulative, PMessageBundle, PSizeInfoMessages;
    vector<MPI_Request> SRequest;  // Requests that are used for each send message
    vector<MPI_Status> SStatus;    // Status of sent messages, used in MPI_Wait
    MilanLongInt MessageIndex = 0; // Pointer for current message
    MilanInt BufferSize;
    MilanLongInt *Buffer;

    // Declare the locks
    omp_lock_t MateLock[NLVer];

    initialize(NLVer, NLEdge, StartIndex,
               EndIndex, &numGhostEdges,
               &numGhostVertices, &S,
               verLocInd, verLocPtr,
               MateLock,
               Ghost2LocalMap, Counter,
               verGhostPtr, verGhostInd,
               tempCounter, GMate,
               Message, QLocalVtx,
               QGhostVtx, QMsgType, QOwner,
               candidateMate, U,
               privateU,
               privateQLocalVtx,
               privateQGhostVtx,
               privateQMsgType,
               privateQOwner);

    finishTime = MPI_Wtime();
    *ph0_time = finishTime - startTime; // Time taken for Phase-0: Initialization

    startTime = MPI_Wtime();

    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////// INITIALIZATION /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // Compute the Initial Matching Set:

    /*
     * OMP PARALLEL_COMPUTE_CANDIDATE_MATE_B has been splitted from
     * PARALLEL_PROCESS_EXPOSED_VERTEX_B in order to better parallelize
     * the two.
     * PARALLEL_COMPUTE_CANDIDATE_MATE_B is now totally parallel.
     */

    PARALLEL_COMPUTE_CANDIDATE_MATE_B(NLVer,
                                      verLocPtr,
                                      verLocInd,
                                      myRank,
                                      edgeLocWeight,
                                      candidateMate);

    /*
     * PARALLEL_PROCESS_EXPOSED_VERTEX_B
     * TODO: write comment
     *
     * TODO: Test when it's actually more efficient to execute this code
     *       in parallel.
     */

    PARALLEL_PROCESS_EXPOSED_VERTEX_B(NLVer,
                                      candidateMate,
                                      verLocInd,
                                      verLocPtr,
                                      StartIndex,
                                      EndIndex,
                                      Mate,
                                      GMate,
                                      Ghost2LocalMap,
                                      edgeLocWeight,
                                      &myCard,
                                      &msgInd,
                                      &NumMessagesBundled,
                                      &S,
                                      verDistance,
                                      PCounter,
                                      Counter,
                                      myRank,
                                      numProcs,
                                      U,
                                      privateU,
                                      QLocalVtx,
                                      QGhostVtx,
                                      QMsgType,
                                      QOwner,
                                      privateQLocalVtx,
                                      privateQGhostVtx,
                                      privateQMsgType,
                                      privateQOwner);

    tempCounter.clear(); // Do not need this any more

    ///////////////////////////////////////////////////////////////////////////////////
    /////////////////////////// PROCESS MATCHED VERTICES //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////

    vector<MilanLongInt> UChunkBeingProcessed;
    UChunkBeingProcessed.reserve(UCHUNK);

    processMatchedVertices(NLVer,
                           UChunkBeingProcessed,
                           U,
                           privateU,
                           StartIndex,
                           EndIndex,
                           &myCard,
                           &msgInd,
                           &NumMessagesBundled,
                           &S,
                           verLocPtr,
                           verLocInd,
                           verDistance,
                           PCounter,
                           Counter,
                           myRank,
                           numProcs,
                           candidateMate,
                           GMate,
                           Mate,
                           Ghost2LocalMap,
                           edgeLocWeight,
                           QLocalVtx,
                           QGhostVtx,
                           QMsgType,
                           QOwner,
                           privateQLocalVtx,
                           privateQGhostVtx,
                           privateQMsgType,
                           privateQOwner);

    /////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////// SEND BUNDLED MESSAGES /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    sendBundledMessages(&numGhostEdges,
                        &BufferSize,
                        Buffer,
                        PCumulative,
                        PMessageBundle,
                        PSizeInfoMessages,
                        PCounter,
                        NumMessagesBundled,
                        &msgActual,
                        &MessageIndex,
                        numProcs,
                        myRank,
                        comm,
                        QLocalVtx,
                        QGhostVtx,
                        QMsgType,
                        QOwner,
                        SRequest,
                        SStatus);

    ///////////////////////// END OF SEND BUNDLED MESSAGES //////////////////////////////////

    finishTime = MPI_Wtime();
    *ph1_time = finishTime - startTime; // Time taken for Phase-1
    *ph1_card = myCard;                 // Cardinality at the end of Phase-1
    startTime = MPI_Wtime();
    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////// MAIN LOOP //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // Main While Loop:
#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << "=========================************===============================" << endl;
    fflush(stdout);
    fflush(stdout);
#endif
#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << ")Entering While(true) loop..";
    fflush(stdout);
#endif
#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << "=========================************===============================" << endl;
    fflush(stdout);
    fflush(stdout);
#endif

    while (true)
    {
#ifdef DEBUG_HANG_
        if (myRank == 0)
            cout << "\n(" << myRank << ") Main loop" << endl;
        fflush(stdout);
#endif
        ///////////////////////////////////////////////////////////////////////////////////
        /////////////////////////// PROCESS MATCHED VERTICES //////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////
        while (/*!Q.empty()*/ !U.empty())
        {
            // Q.pop_front();
            u = U.pop_front(); // Get an element from the queue
#ifdef PRINT_DEBUG_INFO_
            cout << "\n(" << myRank << ")u: " << u;
            fflush(stdout);
#endif
            if ((u >= StartIndex) && (u <= EndIndex))
            { // Process Only If a Local Vertex
                // Get the Adjacency list for u
                adj1 = verLocPtr[u - StartIndex]; // Pointer
                adj2 = verLocPtr[u - StartIndex + 1];
                for (k = adj1; k < adj2; k++)
                {
                    v = verLocInd[k];
                    if ((v >= StartIndex) && (v <= EndIndex))
                    {                                  // v is a Local Vertex:
                        if (Mate[v - StartIndex] >= 0) // v is already matched
                            continue;
#ifdef PRINT_DEBUG_INFO_
                        cout << "\n(" << myRank << ")v: " << v << " c(v)= " << candidateMate[v - StartIndex] << " Mate[v]: " << Mate[v];
                        fflush(stdout);
#endif
                        if (candidateMate[v - StartIndex] == u)
                        { // Only if pointing to the matched vertex
                            // Start: PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
                            // Start: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
                            adj11 = verLocPtr[v - StartIndex];
                            adj12 = verLocPtr[v - StartIndex + 1];
                            w = -1;
                            heaviestEdgeWt = MilanRealMin; // Assign the smallest Value possible first LDBL_MIN
                            for (k1 = adj11; k1 < adj12; k1++)
                            {
                                if ((verLocInd[k1] < StartIndex) || (verLocInd[k1] > EndIndex))
                                {                                                  // Is it a ghost vertex?
                                    if (GMate[Ghost2LocalMap[verLocInd[k1]]] >= 0) // Already matched
                                        continue;
                                }
                                else
                                {                                              // A local vertex
                                    if (Mate[verLocInd[k1] - StartIndex] >= 0) // Already matched
                                        continue;
                                }

                                if ((edgeLocWeight[k1] > heaviestEdgeWt) ||
                                    ((edgeLocWeight[k1] == heaviestEdgeWt) && (w < verLocInd[k1])))
                                {
                                    heaviestEdgeWt = edgeLocWeight[k1];
                                    w = verLocInd[k1];
                                }
                            } // End of for loop
                            candidateMate[v - StartIndex] = w;
                            // End: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
#ifdef PRINT_DEBUG_INFO_
                            cout << "\n(" << myRank << ")" << v << " Points to: " << w;
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
                                    cout << "\n(" << myRank << ")Sending a request message:";
                                    cout << "\n(" << myRank << ")Ghost is " << w << " Owner is: " << findOwnerOfGhost(w, verDistance, myRank, numProcs);
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
                                        // Q.push_back(u);
                                        U.push_back(v);
                                        U.push_back(w);
                                        myCard++;
#ifdef PRINT_DEBUG_INFO_
                                        cout << "\n(" << myRank << ")MATCH: (" << v << "," << w << ") ";
                                        fflush(stdout);
#endif

                                        PROCESS_CROSS_EDGE(Counter, Ghost2LocalMap[w], &S);

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
                                        cout << "\n(" << myRank << ")MATCH: (" << v << "," << w << ") ";
                                        fflush(stdout);
#endif
                                    } // End of if(CandidateMate(w) = v
                                }     // End of Else
                            }         // End of if(w >=0)
                            else
                            { // no dominating edge found: w == -1
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
                                        cout << "\n(" << myRank << ")Ghost is " << w << " Owner is: " << findOwnerOfGhost(w, verDistance, myRank, numProcs);
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
                        } // End of If (candidateMate[v-StartIndex] == u)
                    }     // End of if ( (v >= StartIndex) && (v <= EndIndex) ) //If Local Vertex:
                    else
                    { // Neighbor v is a ghost vertex
                        if (candidateMate[NLVer + Ghost2LocalMap[v]] == u)
                            candidateMate[NLVer + Ghost2LocalMap[v]] = -1;
                        if (v != Mate[u - StartIndex])
                        { // u is a local vertex
                            // Build the Message Packet:
                            Message[0] = u;       // LOCAL
                            Message[1] = v;       // GHOST
                            Message[2] = SUCCESS; // TYPE
                            // Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
                            cout << "\n(" << myRank << ")Sending a success message: ";
                            cout << "\n(" << myRank << ")Ghost is " << v << " Owner is: " << findOwnerOfGhost(v, verDistance, myRank, numProcs);
                            fflush(stdout);
#endif
                            ghostOwner = findOwnerOfGhost(v, verDistance, myRank, numProcs);
                            assert(ghostOwner != -1);
                            assert(ghostOwner != myRank);
                            MPI_Bsend(&Message[0], 3, TypeMap<MilanLongInt>(), ghostOwner, ComputeTag, comm);
                            msgInd++;
                            msgActual++;
#ifdef DEBUG_GHOST_
                            if ((u < StartIndex) || (u > EndIndex))
                            {
                                cout << "\n(" << myRank << ") " << __LINE__ << " From Send: should not happen: u= " << u << " v= " << v << " StartIndex " << StartIndex << " EndIndex " << EndIndex << endl;
                                fflush(stdout);
                            }
#endif

                        } // End of If( v != Mate[u] )
                    }     // End of Else //A Ghost Vertex
                }         // End of For Loop adj(u)
            }             // End of if ( (u >= StartIndex) && (u <= EndIndex) ) //Process Only If a Local Vertex
        }                 // End of while ( /*!Q.empty()*/ !U.empty() )
        ///////////////////////// END OF PROCESS MATCHED VERTICES /////////////////////////

        //// BREAK IF NO MESSAGES EXPECTED /////////
#ifdef PRINT_DEBUG_INFO_
        cout << "\n(" << myRank << ")Deciding whether to break: S= " << S << endl;
#endif

        if (S == 0)
        {
#ifdef DEBUG_HANG_
            cout << "\n(" << myRank << ") Breaking out" << endl;
            fflush(stdout);
#endif
            break;
        }
        ///////////////////////////////////////////////////////////////////////////////////
        /////////////////////////// PROCESS MESSAGES //////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////

        processMessages(NLVer,
                        Mate,
                        candidateMate,
                        Ghost2LocalMap,
                        GMate,
                        Counter,
                        StartIndex,
                        EndIndex,
                        &myCard,
                        &msgInd,
                        &msgActual,
                        edgeLocWeight,
                        verDistance,
                        verLocPtr,
                        k,
                        verLocInd,
                        numProcs,
                        myRank,
                        comm,
                        Message,
                        numGhostEdges,
                        u,
                        v,
                        &S,
                        U);

        ///////////////////////// END OF PROCESS MESSAGES /////////////////////////////////
#ifdef PRINT_DEBUG_INFO_
        cout << "\n(" << myRank << ")Finished Message processing phase: S= " << S;
        fflush(stdout);
        cout << "\n(" << myRank << ")** SENT     : ACTUAL= " << msgActual;
        fflush(stdout);
        cout << "\n(" << myRank << ")** SENT     : INDIVIDUAL= " << msgInd << endl;
        fflush(stdout);
#endif
    } // End of while (true)

    clean(NLVer,
          myRank,
          MessageIndex,
          SRequest,
          SStatus,
          BufferSize,
          Buffer,
          msgActual,
          msgActualSent,
          msgInd,
          msgIndSent,
          NumMessagesBundled,
          msgPercent,
          MateLock);

    finishTime = MPI_Wtime();
    *ph2_time = finishTime - startTime; // Time taken for Phase-2
    *ph2_card = myCard;                 // Cardinality at the end of Phase-2
}
// End of algoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMate
#endif

#endif