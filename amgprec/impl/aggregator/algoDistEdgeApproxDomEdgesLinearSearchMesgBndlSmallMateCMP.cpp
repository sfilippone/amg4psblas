#include "MatchBoxPC.h"
#include <omp.h>
#include <stdio.h>
#include "isAlreadyMatched.cpp"
#include "findOwnerOfGhost.cpp"
#include "computeCandidateMate.cpp"
#include "initialize.cpp"

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

#define UCHUNK 1000

#ifdef SERIAL_MPI
#else
//MPI type map
template<typename T> MPI_Datatype TypeMap();
template<> inline MPI_Datatype TypeMap<int64_t>() { return MPI_LONG_LONG; }
template<> inline MPI_Datatype TypeMap<int>() { return MPI_INT; }
template<> inline MPI_Datatype TypeMap<double>() { return MPI_DOUBLE; }
template<> inline MPI_Datatype TypeMap<float>() { return MPI_FLOAT; }

// DOUBLE PRECISION VERSION
//WARNING: The vertex block on a given rank is contiguous
void dalgoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMateCMP(
        MilanLongInt NLVer, MilanLongInt NLEdge,
        MilanLongInt* verLocPtr, MilanLongInt* verLocInd,
        MilanReal* edgeLocWeight,
        MilanLongInt* verDistance,
        MilanLongInt* Mate,
        MilanInt myRank, MilanInt numProcs, MPI_Comm comm,
        MilanLongInt* msgIndSent, MilanLongInt* msgActualSent,
        MilanReal* msgPercent,
        MilanReal* ph0_time, MilanReal* ph1_time, MilanReal* ph2_time,
        MilanLongInt* ph1_card, MilanLongInt* ph2_card ) {

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
    cout<<"\n("<<myRank<<")Within algoEdgeApproxDominatingEdgesLinearSearchMessageBundling()"; fflush(stdout);
#endif

#ifdef PRINT_DEBUG_INFO_
    cout<<"\n("<<myRank<<") verDistance ["<< verDistance[0] << "," << verDistance[1] << "," << verDistance[2] <<"," << verDistance[3] <<"]"; fflush(stdout);
#endif
#ifdef DEBUG_HANG_
    if (myRank == 0) cout<<"\n("<<myRank<<") verDistance ["<< verDistance[0] << "," << verDistance[1] << "," << verDistance[2] <<"," << verDistance[3] <<"]"; fflush(stdout);
#endif

    //inputSubGraph.getStartEndIndices(StartIndex, EndIndex);
    MilanLongInt StartIndex = verDistance[myRank]; //The starting vertex owned by the current rank
    //MilanLongInt EndIndex = verDistance[myRank+1]; //The ending vertex owned by the current rank
    MilanLongInt EndIndex = verDistance[myRank + 1] - 1; //The ending vertex owned by the current rank

    MPI_Status computeStatus;
    const int ComputeTag = 7;  //Predefined tag
    const int BundleTag = 9;   //Predefined tag
    int error_codeC;
    error_codeC = MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    char error_message[MPI_MAX_ERROR_STRING];
    int message_length;

    //MilanLongInt NLVer=0, NLEdge=0, StartIndex=0, EndIndex=0;
    MilanLongInt msgActual = 0, msgInd = 0;
    MilanReal heaviestEdgeWt = 0.0f; //Assumes positive weight
    MilanReal startTime, finishTime;
    //MilanReal Precision = MPI_Wtick(); //Get the precision of the MPI Timer
    startTime = MPI_Wtime();

    //Data structures for sending and receiving messages:
    vector<MilanLongInt> Message; // [ u, v, message_type ]
    Message.resize(3,-1);
    const MilanLongInt REQUEST  = 1;
    const MilanLongInt SUCCESS  = 2;
    const MilanLongInt FAILURE  = 3;
    const MilanLongInt SIZEINFO = 4;
    MilanLongInt message_type = 0;
    //Data structures for Message Bundling:
    //Although up to two messages can be sent along any cross edge,
    //only one message will be sent in the initialization phase -
    //one of: REQUEST/FAILURE/SUCCESS
    vector<MilanLongInt> QLocalVtx, QGhostVtx, QMsgType;
    vector<MilanInt> QOwner; // Changed by Fabio to be an integer, addresses needs to be integers!

    MilanLongInt* PCounter = new MilanLongInt [numProcs];
    for (int i = 0; i < numProcs; i++)
        PCounter[i] = 0;


    MilanLongInt NumMessagesBundled = 0;
    MilanInt ghostOwner = 0; // Changed by Fabio to be an integer, addresses needs to be integers!
    MilanLongInt* candidateMate = nullptr;
#ifdef PRINT_DEBUG_INFO_
    cout<<"\n("<<myRank<<")NV: "<<NLVer<<"  Edges: "<<NLEdge; fflush(stdout);
    cout<<"\n("<<myRank<<")StartIndex: "<<StartIndex<<"  EndIndex: "<<EndIndex; fflush(stdout);
#endif
    //Other Variables:
    MilanLongInt u = -1, v = -1, w = -1, i = 0;
    MilanLongInt k = -1, adj1 = -1, adj2 = -1;
    MilanLongInt k1 = -1, adj11 = -1, adj12 = -1;
    MilanLongInt myCard = 0;
    MilanInt Sender = 0; // This is the rank of the sending nodes, it has to be an integer! Fabio

    //Build the Ghost Vertex Set: Vg
    map <MilanLongInt, MilanLongInt> Ghost2LocalMap; //Map each ghost vertex to a local vertex
    vector <MilanLongInt> Counter;  //Store the edge count for each ghost vertex
    MilanLongInt numGhostVertices = 0, numGhostEdges = 0; //Number of Ghost vertices

#ifdef PRINT_DEBUG_INFO_
    cout<<"\n("<<myRank<<")About to compute Ghost Vertices..."; fflush(stdout);
#endif
#ifdef DEBUG_HANG_
    if (myRank == 0)     cout<<"\n("<<myRank<<")About to compute Ghost Vertices..."; fflush(stdout);
#endif

    //Define Adjacency Lists for Ghost Vertices:
    //cout<<"Building Ghost data structures ... \n\n";
    vector <MilanLongInt> verGhostPtr, verGhostInd, tempCounter;
    //Mate array for ghost vertices:
    vector <MilanLongInt> GMate;  //Proportional to the number of ghost vertices
    MilanLongInt S;
    MilanLongInt privateMyCard = 0;
    staticQueue U, privateU, privateQLocalVtx, privateQGhostVtx, privateQMsgType, privateQOwner;
    MilanLongInt myIndex = 0;
    vector <MilanLongInt> PCumulative, PMessageBundle, PSizeInfoMessages;
    vector <MPI_Request> SRequest; //Requests that are used for each send message
    vector <MPI_Status> SStatus;   //Status of sent messages, used in MPI_Wait
    MilanLongInt MessageIndex = 0; //Pointer for current message
    MilanInt OneMessageSize = 0;
    MilanLongInt numMessagesToSend;
    MilanInt BufferSize;
    MilanLongInt *Buffer;
    bool isEmpty;

    //Declare the locks
    // TODO destroy the locks
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
                privateQOwner
                );
                        
    finishTime = MPI_Wtime();
    *ph0_time = finishTime - startTime; //Time taken for Phase-0: Initialization      

                    
    startTime = MPI_Wtime();
    

    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////// INITIALIZATION /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    //Compute the Initial Matching Set:

#pragma omp parallel private(k, u, w, v, k1, adj1, adj2, adj11, adj12, heaviestEdgeWt, ghostOwner, privateMyCard, isEmpty) firstprivate(privateU, StartIndex, EndIndex, privateQLocalVtx, privateQGhostVtx, privateQMsgType, privateQOwner) default(shared) num_threads(4)
    {
        /*
        * OMP PARALLEL_COMPUTE_CANDIDATE_MATE_B has been splitted from
        * PARALLEL_PROCESS_EXPOSED_VERTEX_B in order to better parallelize
        * the two.
        * In particular PARALLEL_COMPUTE_CANDIDATE_MATE_B is now totally parallel.
        */

#pragma omp for schedule(static)
        for ( v=0; v < NLVer; v++ ) {
#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")Processing: "<<v+StartIndex<<endl; fflush(stdout);
#endif
            //Start: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
            candidateMate[v] = firstComputeCandidateMate(verLocPtr[v], verLocPtr[v + 1], verLocInd, edgeLocWeight);
            //End: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
        }

        /*
         * PARALLEL_PROCESS_EXPOSED_VERTEX_B
         * The sequential version could be a bit more
         * efficient.
         *
         * TODO: Maybe it is possible to append the values of QLocalVtx, QGhostVtx, QMsgType and QOwner
         *       first in a local variable and then, only at the end, append them to the real data structure
         *       to remove the critical sections.
         *
         * TODO: Test when it's more efficient to execute this code
         *       in parallel.
         */


#pragma omp for reduction(+: msgInd, NumMessagesBundled, myCard, PCounter[:numProcs]) schedule(static)
        for (v = 0; v < NLVer; v++) {
            //Start: PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
            k = candidateMate[v];
            candidateMate[v] = verLocInd[k];
            w = candidateMate[v];

#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")Processing: "<<v+StartIndex<<endl; fflush(stdout);
#endif

#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")"<<v+StartIndex<<" Points to: "<<w; fflush(stdout);
#endif
            //If found a dominating edge:
            if (w >= 0) {

                if (isAlreadyMatched(verLocInd[k], StartIndex, EndIndex, GMate, Mate, Ghost2LocalMap)) {
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

                if (w >= 0) {

                    myCard++;
                    if ((w < StartIndex) || (w > EndIndex)) { //w is a ghost vertex
#ifdef PRINT_DEBUG_INFO_
                        cout<<"\n("<<myRank<<")Sending a request message (291):";
                        cout<<"\n("<<myRank<<")Local is: "<<v+StartIndex<<" Ghost is "<<w<<" Owner is: "<< findOwnerOfGhost(w, verDistance, myRank, numProcs) <<endl;
                        fflush(stdout);
#endif

                        msgInd++;
                        NumMessagesBundled++;
                        ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
                        assert(ghostOwner != -1);
                        assert(ghostOwner != myRank);
                        PCounter[ghostOwner]++;

                        /*
                        //TODO why does it fail if I use a private data structure???
                        privateQLocalVtx.push_back(v + StartIndex);
                        privateQGhostVtx.push_back(w);
                        privateQMsgType.push_back(REQUEST);
                        privateQOwner.push_back(ghostOwner);
                        */

#pragma omp critical(MSG)
                        {

                            QLocalVtx.push_back(v + StartIndex);
                            QGhostVtx.push_back(w);
                            QMsgType.push_back(REQUEST);
                            QOwner.push_back(ghostOwner);
                        } // end of critical region
                    

                        if (candidateMate[NLVer + Ghost2LocalMap[w]] == v + StartIndex) {

                            privateU.push_back(v + StartIndex);
                            privateU.push_back(w);
                            Mate[v] = w;
                            //FIXME could this instruction create errors?
                            GMate[Ghost2LocalMap[w]] = v + StartIndex; //w is a Ghost

#ifdef PRINT_DEBUG_INFO_
                            cout<<"\n("<<myRank<<")MATCH: ("<<v+StartIndex<<","<<w<<")"; fflush(stdout);
#endif
                            //Decrement the counter:
                            //Start: PARALLEL_PROCESS_CROSS_EDGE_B(v)
#pragma omp critical
                            {
                                if (Counter[Ghost2LocalMap[w]] > 0) {

                                    Counter[Ghost2LocalMap[w]] -= 1; //Decrement
                                    if (Counter[Ghost2LocalMap[w]] == 0) {
                                        S--; //Decrement S
#ifdef PRINT_DEBUG_INFO_
                                        cout<<"\n("<<myRank<<")Decrementing S: Ghost vertex "<<w<<" has received all its messages";
                                            fflush(stdout);
#endif
                                    }
                                }
                            } //End of if Counter[w] > 0
                            //End: PARALLEL_PROCESS_CROSS_EDGE_B(v)
                        } //End of if CandidateMate[w] = v


                    } //End of if a Ghost Vertex
                    else { // w is a local vertex

                        if (candidateMate[w - StartIndex] == (v + StartIndex)) {
                            privateU.push_back(v + StartIndex);
                            privateU.push_back(w);

                            Mate[v] = w;  //v is local
                            //FIXME this instruction could create errors
                            Mate[w - StartIndex] = v + StartIndex; //w is local


#ifdef PRINT_DEBUG_INFO_
                            cout<<"\n("<<myRank<<")MATCH: ("<<v+StartIndex<<","<<w<<") "; fflush(stdout);
#endif

                        } //End of if ( candidateMate[w-StartIndex] == (v+StartIndex) )
                    } //End of Else

                    continue;
                } //End of second if

            } //End of if(w >=0)

            //This piece of code is executed a really small amount of times, I will not allocate a
            //huge amount of memory to the private data structures.
            adj11 = verLocPtr[v];
            adj12 = verLocPtr[v + 1];
            for (k1 = adj11; k1 < adj12; k1++) {
                w = verLocInd[k1];
                if ((w < StartIndex) || (w > EndIndex)) { //A ghost

#ifdef PRINT_DEBUG_INFO_
                    cout<<"\n("<<myRank<<")Sending a failure message: ";
                cout<<"\n("<<myRank<<")Ghost is "<<w<<" Owner is: "<<findOwnerOfGhost(w, verDistance, myRank, numProcs);
                fflush(stdout);
#endif

                    msgInd++;
                    NumMessagesBundled++;
                    ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
                    assert(ghostOwner != -1);
                    assert(ghostOwner != myRank);
                    PCounter[ghostOwner]++;
                    QLocalVtx.push_back(v + StartIndex);
                    QGhostVtx.push_back(w);
                    QMsgType.push_back(FAILURE);
                    QOwner.push_back(ghostOwner);

                } //End of if(GHOST)
            } //End of for loop
            //End:   PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
        } //End of for ( v=0; v < NLVer; v++ )


        #pragma omp critical(privateMsg)
        {
            while (!privateQLocalVtx.empty()) {
                
                QLocalVtx.push_back(privateQLocalVtx.pop_front());
                QGhostVtx.push_back(privateQGhostVtx.pop_front());
                QMsgType.push_back(privateQMsgType.pop_front());
                QOwner.push_back(privateQOwner.pop_front());

            }

        }

#pragma omp critical(U)
        {
            while (!privateU.empty())
            {
                U.push_back(privateU.pop_front());
            }
        }

#pragma omp master
        {
            tempCounter.clear(); //Do not need this any more
        }

#ifdef PRINT_DEBUG_INFO_
        cout<<"\n("<<myRank<<"=========================************==============================="<<endl; fflush(stdout);
    fflush(stdout);
#endif
        ///////////////////////////////////////////////////////////////////////////////////
        /////////////////////////// PROCESS MATCHED VERTICES //////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////
        isEmpty = false;

#ifdef COUNT_LOCAL_VERTEX
        MilanLongInt localVertices = 0;
#endif

        //TODO what would be the optimal UCHUNK
        vector <MilanLongInt> Us;
        Us.reserve(UCHUNK);

        while( true ) {

            Us.clear();
#pragma omp critical(U)
            {
                //If U is emptu and there are no new node to add to U
                if (U.empty() && privateU.empty())
                    isEmpty = true;
                else {
                    if (U.empty() && !privateU.empty()) // If U is empty but there are nodes in private U
                        while (!privateU.empty()) {
                            U.push_back(privateU.pop_front());
                            myCard += privateMyCard;
                        }
                    for (int i = 0; i < UCHUNK; i++) { // Pop the new nodes
                        if (U.empty()) break;
                        Us.push_back(U.pop_front());
                    }
                }
            } // End of critical U
            if (isEmpty) break;

            for (MilanLongInt u : Us)
            {
#ifdef PRINT_DEBUG_INFO_
                cout<<"\n("<<myRank<<")u: "<<u; fflush(stdout);
#endif
                if ((u >= StartIndex) && (u <= EndIndex)) { //Process Only the Local Vertices

#ifdef COUNT_LOCAL_VERTEX
                    localVertices ++;
#endif

                    //Get the Adjacency list for u
                    adj1 = verLocPtr[u - StartIndex];  //Pointer
                    adj2 = verLocPtr[u - StartIndex + 1];
                    for (k = adj1; k < adj2; k++) {
                        v = verLocInd[k];

                        if ((v >= StartIndex) && (v <= EndIndex)) { //If Local Vertex:
#pragma omp critical(innerProcessMatched)
                            {

#ifdef PRINT_DEBUG_INFO_
                                cout<<"\n("<<myRank<<")v: "<<v<<" c(v)= "<<candidateMate[v-StartIndex]<<" Mate[v]: "<<Mate[v];
                        fflush(stdout);
#endif


                                //If the current vertex is pointing to a matched vertex and is not matched
                                //FIXME is there a way to make candidateMate private?
                                //      for the moment it could generate an error.
                                if (not isAlreadyMatched(v, StartIndex, EndIndex, GMate, Mate, Ghost2LocalMap) and
                                    candidateMate[v - StartIndex] == u) {


                                    //Start: PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
                                    //Start: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
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

                                    //End: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
#ifdef PRINT_DEBUG_INFO_
                                    cout<<"\n("<<myRank<<")"<<v<<" Points to: "<<w; fflush(stdout);
#endif
                                    //If found a dominating edge:
                                    if (w >= 0) {

                                        //TODO is it possible to lock without a critical region?
                                        //TODO there must be a more elegant and efficient way to do this
                                        /*
                                        while(true) {
                                            if (omp_test_lock(&MateLock[v - StartIndex])) {
                                                if (omp_test_lock(&MateLock[w - StartIndex])) break;
                                                else omp_unset_lock(&MateLock[v - StartIndex]);
                                            }
                                        }
                                        */


                                        if ((w < StartIndex) || (w > EndIndex)) { //A ghost
#ifdef PRINT_DEBUG_INFO_
                                            cout<<"\n("<<myRank<<")Sending a request message:";
                                    cout<<"\n("<<myRank<<")Ghost is "<<w<<" Owner is: "<<findOwnerOfGhost(w, verDistance, myRank, numProcs);
#endif

                                            QLocalVtx.push_back(v);
                                            QGhostVtx.push_back(w);
                                            QMsgType.push_back(REQUEST);
                                            ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
                                            assert(ghostOwner != -1);
                                            assert(ghostOwner != myRank);
                                            QOwner.push_back(ghostOwner);
                                            PCounter[ghostOwner]++;
                                            NumMessagesBundled++;
                                            msgInd++;
                                            if (candidateMate[NLVer + Ghost2LocalMap[w]] == v) {
                                                Mate[v - StartIndex] = w;  //v is a local vertex
                                                GMate[Ghost2LocalMap[w]] = v;  //w is a ghost vertex
                                                //Q.push_back(u);
                                                privateU.push_back(v);
                                                privateU.push_back(w);
                                                privateMyCard++;
#ifdef PRINT_DEBUG_INFO_
                                                cout<<"\n("<<myRank<<")MATCH: ("<<v<<","<<w<<") "; fflush(stdout);
#endif
                                                //Decrement the counter:
                                                //Start: PARALLEL_PROCESS_CROSS_EDGE_B(v,w)
                                                if (Counter[Ghost2LocalMap[w]] > 0) {
                                                    Counter[Ghost2LocalMap[w]] = Counter[Ghost2LocalMap[w]] - 1; //Decrement
                                                    if (Counter[Ghost2LocalMap[w]] == 0) {
                                                        S--; //Decrement S
#ifdef PRINT_DEBUG_INFO_
                                                        cout<<"\n("<<myRank<<")Decrementing S: Ghost vertex "<<w<<" has received all its messages";
                                                fflush(stdout);
#endif
                                                    }
                                                } //End of if Counter[w] > 0
                                                //End: PARALLEL_PROCESS_CROSS_EDGE_B(v,w)
                                            } //End of if CandidateMate[w] = v
                                        } //End of if a Ghost Vertex
                                        else { //w is a local vertex
                                            if (candidateMate[w - StartIndex] == v) {
                                                Mate[v - StartIndex] = w;  //v is a local vertex
                                                Mate[w - StartIndex] = v;  //w is a local vertex
                                                //Q.push_back(u);
                                                privateU.push_back(v);
                                                privateU.push_back(w);
                                                privateMyCard++;
#ifdef PRINT_DEBUG_INFO_
                                                cout<<"\n("<<myRank<<")MATCH: ("<<v<<","<<w<<") "; fflush(stdout);
#endif
                                            } //End of if(CandidateMate(w) = v
                                        } //End of Else

                                        //omp_unset_lock(&MateLock[v - StartIndex]);
                                        //omp_unset_lock(&MateLock[w - StartIndex]);

                                    } //End of if(w >=0)
                                    else {
                                        adj11 = verLocPtr[v - StartIndex];
                                        adj12 = verLocPtr[v - StartIndex + 1];
                                        for (k1 = adj11; k1 < adj12; k1++) {
                                            w = verLocInd[k1];
                                            if ((w < StartIndex) || (w > EndIndex)) { //A ghost

#ifdef PRINT_DEBUG_INFO_
                                                cout<<"\n("<<myRank<<")Sending a failure message: ";
                                        cout<<"\n("<<myRank<<")Ghost is "<<w<<" Owner is: "<<findOwnerOfGhost(w, verDistance, myRank, numProcs);
                                        fflush(stdout);
#endif
                                                /* MPI_Bsend(&Message[0], 3, MPI_INT, inputSubGraph.findOwner(w),
                                                 ComputeTag, comm); */
                                                QLocalVtx.push_back(v);
                                                QGhostVtx.push_back(w);
                                                QMsgType.push_back(FAILURE);
                                                //ghostOwner = inputSubGraph.findOwner(w);
                                                ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
                                                assert(ghostOwner != -1);
                                                assert(ghostOwner != myRank);
                                                QOwner.push_back(ghostOwner);
                                                PCounter[ghostOwner]++;
                                                NumMessagesBundled++;
                                                msgInd++;
                                            } //End of if(GHOST)
                                        } //End of for loop
                                    } // End of Else: w == -1
                                    //End:   PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)

                                } //End of If (candidateMate[v-StartIndex] == u

                            } //End of critical region if

                        } //End of if ( (v >= StartIndex) && (v <= EndIndex) ) //If Local Vertex:
                        else { //Neighbor is a ghost vertex

#pragma omp critical(innerProcessMatched)
                            {

                                //while(!omp_test_lock(&MateLock[u - StartIndex]));

                                if (candidateMate[NLVer + Ghost2LocalMap[v]] == u)
                                    candidateMate[NLVer + Ghost2LocalMap[v]] = -1;
                                if (v != Mate[u - StartIndex]) { //u is local
                                    //Build the Message Packet:
                                    //Message[0] = u; //LOCAL
                                    //Message[1] = v; //GHOST
                                    //Message[2] = SUCCESS;  //TYPE
                                    //Send a Request (Asynchronous)

#ifdef PRINT_DEBUG_INFO_
                                    cout<<"\n("<<myRank<<")Sending a success message: ";
                            cout<<"\n("<<myRank<<")Ghost is "<<v<<" Owner is: "<<findOwnerOfGhost(v, verDistance, myRank, numProcs)<<"\n"; fflush(stdout);
#endif

                                    QLocalVtx.push_back(u);
                                    QGhostVtx.push_back(v);
                                    QMsgType.push_back(SUCCESS);
                                    ghostOwner = findOwnerOfGhost(v, verDistance, myRank, numProcs);
                                    assert(ghostOwner != -1);
                                    assert(ghostOwner != myRank);
                                    QOwner.push_back(ghostOwner);
                                    PCounter[ghostOwner]++;
                                    NumMessagesBundled++;
                                    msgInd++;
                                } //End of If( v != Mate[u] )

                                //omp_unset_lock(&MateLock[u - StartIndex]);

                            } //End of critical region
                        } //End of Else //A Ghost Vertex

                    } //End of For Loop adj(u)

                } //End of if ( (u >= StartIndex) && (u <= EndIndex) ) //Process Only If a Local Vertex

                //Avoid to ask for the critical section if there is nothing to add
                if (privateU.size() < UCHUNK && !U.empty()) continue;
#pragma omp critical(U)
                {
                    while (!privateU.empty()) {
                        U.push_back(privateU.pop_front());
                    }

                    myCard += privateMyCard;
                } //End of critical U

            }
        } //End of while ( /*!Q.empty()*/ !U.empty() )

        #pragma omp critical(privateMsg)
        {
            while (!privateQLocalVtx.empty()) {

                QLocalVtx.push_back(privateQLocalVtx.pop_front());
                QGhostVtx.push_back(privateQGhostVtx.pop_front());
                QMsgType.push_back(privateQMsgType.pop_front());
                QOwner.push_back(privateQOwner.pop_front());

            }

        }


#ifdef COUNT_LOCAL_VERTEX
        printf("Count local vertexes: %ld for thread %d of processor %d\n",
              localVertices,
              omp_get_thread_num(),
              myRank);
#endif


        ///////////////////////// END OF PROCESS MATCHED VERTICES /////////////////////////
#ifdef DEBUG_HANG_
        if (myRank == 0) cout<<"\n("<<myRank<<") Send Bundles" <<endl; fflush(stdout);
#endif
        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////// SEND BUNDLED MESSAGES /////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////
#pragma omp barrier
#pragma omp master
        {
            //Data structures for Bundled Messages:
            try {
                PMessageBundle.reserve(NumMessagesBundled * 3); //Three integers per message
                PCumulative.reserve(numProcs + 1); //Similar to Row Pointer vector in CSR data structure
                PSizeInfoMessages.reserve(numProcs * 3); //Buffer to hold the Size info message packets
            } catch (length_error) {
                cout << "Error in function algoDistEdgeApproxDominatingEdgesMessageBundling: \n";
                cout << "Not enough memory to allocate the internal variables \n";
                exit(1);
            }
            PMessageBundle.resize(NumMessagesBundled * 3, -1);//Initialize
            PCumulative.resize(numProcs + 1, 0); //Only initialize the counter variable
            PSizeInfoMessages.resize(numProcs * 3, 0);


            for (MilanInt i = 0; i < numProcs; i++) // Changed by Fabio to be an integer, addresses needs to be integers!
                PCumulative[i + 1] = PCumulative[i] + PCounter[i];

            //OMP not worth parallelizing
            //Reuse PCounter to keep track of how many messages were inserted:
            for (MilanInt i = 0; i < numProcs; i++) // Changed by Fabio to be an integer, addresses needs to be integers!
                PCounter[i] = 0;
            //Build the Message Bundle packet:

            //OMP Not parallelizable
            for (MilanInt i=0; i<NumMessagesBundled; i++) { // Changed by Fabio to be an integer, addresses needs to be integers!
                myIndex = ( PCumulative[QOwner[i]] + PCounter[QOwner[i]] )*3;
                PMessageBundle[myIndex+0] = QLocalVtx[i];
                PMessageBundle[myIndex+1] = QGhostVtx[i];
                PMessageBundle[myIndex+2] = QMsgType[i];
                PCounter[QOwner[i]]++;
            }

            //Send the Bundled Messages: Use ISend

            try {
                SRequest.reserve(numProcs * 2); //At most two messages per processor
                SStatus.reserve(numProcs * 2);//At most two messages per processor
            } catch (length_error) {
                cout << "Error in function algoDistEdgeApproxDominatingEdgesLinearSearchImmediateSend: \n";
                cout << "Not enough memory to allocate the internal variables \n";
                exit(1);
            }
            MPI_Request myReq; //A sample request
            SRequest.resize(numProcs * 2, myReq);
            MPI_Status myStat; //A sample status
            SStatus.resize(numProcs * 2, myStat);

            //Send the Messages
            for (MilanInt i = 0; i < numProcs; i++) { // Changed by Fabio to be an integer, addresses needs to be integers!
                if (i == myRank) //Do not send anything to yourself
                    continue;
                //Send the Message with information about the size of next message:
                //Build the Message Packet:
                PSizeInfoMessages[i * 3 + 0] = (PCumulative[i + 1] - PCumulative[i]) * 3; // # of integers in the next message
                PSizeInfoMessages[i * 3 + 1] = -1; //Dummy packet
                PSizeInfoMessages[i * 3 + 2] = SIZEINFO;  //TYPE
                //Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
                cout<<"\n("<<myRank<<")Sending bundled message to process "<<i<<" size: "<<PSizeInfoMessages[i*3+0]<<endl;
                fflush(stdout);
#endif
                if (PSizeInfoMessages[i * 3 + 0] > 0) { //Send only if it is a nonempty packet
                    MPI_Isend(&PSizeInfoMessages[i * 3 + 0], 3, TypeMap<MilanLongInt>(), i, ComputeTag, comm,
                              &SRequest[MessageIndex]);
                    msgActual++;
                    MessageIndex++;
                    //Now Send the message with the data packet:
#ifdef PRINT_DEBUG_INFO_
                    cout<<"\n("<<myRank<<")Sending Bundle to : "<<i<<endl;
                    for (k=(PCumulative[i]*3); k< (PCumulative[i]*3+PSizeInfoMessages[i*3+0]); k++)
                        cout<<PMessageBundle[k]<<",";
                    cout<<endl;
                    fflush(stdout);
#endif
                    MPI_Isend(&PMessageBundle[PCumulative[i] * 3], PSizeInfoMessages[i * 3 + 0],
                              TypeMap<MilanLongInt>(), i, BundleTag, comm, &SRequest[MessageIndex]);
                    MessageIndex++;
                } //End of if size > 0
            }
            //Free up temporary memory:
            PCumulative.clear();
            QLocalVtx.clear();
            QGhostVtx.clear();
            QMsgType.clear();
            QOwner.clear();


#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")Number of Ghost edges = "<<numGhostEdges;
    cout<<"\n("<<myRank<<")Total number of potential message X 2 = "<<numGhostEdges*2;
    cout<<"\n("<<myRank<<")Number messages already sent in bundles = "<<NumMessagesBundled;
    if (numGhostEdges>0) {
      cout<<"\n("<<myRank<<")Percentage of total = "<<((double)NumMessagesBundled/(double)(numGhostEdges*2))*100.0<<"% \n";
    }
    fflush(stdout);
#endif

            //Allocate memory for MPI Send messages:
            /* WILL COME BACK HERE - NO NEED TO STORE ALL THIS MEMORY !! */
            OneMessageSize=0;
            MPI_Pack_size(3, TypeMap<MilanLongInt>(), comm, &OneMessageSize); //Size of one message packet
            //How many messages to send?
            //Potentially three kinds of messages will be sent/received:
            //Request, Success, Failure.
            //But only two will be sent from a given processor.
            //Substract the number of messages that have already been sent as bundled messages:
            numMessagesToSend = numGhostEdges*2 - NumMessagesBundled;
            BufferSize = (OneMessageSize+MPI_BSEND_OVERHEAD)*numMessagesToSend;

            Buffer=0;
#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")Size of One Message from PACK= "<<OneMessageSize;
    cout<<"\n("<<myRank<<")Size of Message overhead = "<<MPI_BSEND_OVERHEAD;
    cout<<"\n("<<myRank<<")Number of Ghost edges = "<<numGhostEdges;
    cout<<"\n("<<myRank<<")Number of remaining message = "<<numMessagesToSend;
    cout<<"\n("<<myRank<<")BufferSize = "<<BufferSize;
    cout<<"\n("<<myRank<<")Attaching Buffer on.. ";
    fflush(stdout);
#endif
            if ( BufferSize > 0 ) {
                Buffer = (MilanLongInt *) malloc(BufferSize);  //Allocate memory
                if ( Buffer == 0 ) {
                    cout<<"Error in function algoDistEdgeApproxDominatingEdgesLinearSearch: \n";
                    cout<<"Not enough memory to allocate for send buffer on process "<<myRank<<"\n";
                    exit(1);
                }
                MPI_Buffer_attach(Buffer, BufferSize); //Attach the Buffer
            }
        } //End of master

    } // end of parallel region
    ///////////////////////// END OF SEND BUNDLED MESSAGES //////////////////////////////////

    finishTime = MPI_Wtime();
    *ph1_time = finishTime-startTime; //Time taken for Phase-1
    *ph1_card = myCard; //Cardinality at the end of Phase-1
    startTime = MPI_Wtime();
    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////// MAIN LOOP //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    //Main While Loop:
#ifdef PRINT_DEBUG_INFO_
    cout<<"\n("<<myRank<<"=========================************==============================="<<endl; fflush(stdout);
    fflush(stdout);
#endif
#ifdef PRINT_DEBUG_INFO_
    cout<<"\n("<<myRank<<")Entering While(true) loop.."; fflush(stdout);
    //U.display(); fflush(stdout);
#endif
#ifdef PRINT_DEBUG_INFO_
    cout<<"\n("<<myRank<<"=========================************==============================="<<endl; fflush(stdout);
    fflush(stdout);
#endif
    //Buffer to receive bundled messages
    //Maximum messages that can be received from any processor is
    //twice the edge cut: REQUEST; REQUEST+(FAILURE/SUCCESS)
    vector<MilanLongInt> ReceiveBuffer;
    MilanLongInt bundleSize=0, bundleCounter=0;
    try {
        ReceiveBuffer.reserve(numGhostEdges*2*3); //Three integers per cross edge
    } catch ( length_error ) {
        cout<<"Error in function algoDistEdgeApproxDominatingEdgesMessageBundling: \n";
        cout<<"Not enough memory to allocate the internal variables \n";
        exit(1);
    }
    while ( true ) {
#ifdef DEBUG_HANG_
        if (myRank == 0) cout<<"\n("<<myRank<<") Main loop" <<endl; fflush(stdout);
#endif
        ///////////////////////////////////////////////////////////////////////////////////
        /////////////////////////// PROCESS MATCHED VERTICES //////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////
        while ( /*!Q.empty()*/ !U.empty() ) {
            //Q.pop_front();
            u = U.pop_front(); //Get an element from the queue
#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")u: "<<u; fflush(stdout);
#endif
            if ( (u >= StartIndex) && (u <= EndIndex) ) { //Process Only If a Local Vertex
                //Get the Adjacency list for u
                adj1 = verLocPtr[u-StartIndex];  //Pointer
                adj2 = verLocPtr[u-StartIndex+1];
                for( k = adj1; k < adj2; k++ ) {
                    v = verLocInd[k];
                    if ( (v >= StartIndex) && (v <= EndIndex) ) { //v is a Local Vertex:
                        if ( Mate[v-StartIndex] >= 0 )   // v is already matched
                            continue;
#ifdef PRINT_DEBUG_INFO_
                        cout<<"\n("<<myRank<<")v: "<<v<<" c(v)= "<<candidateMate[v-StartIndex]<<" Mate[v]: "<<Mate[v];
                        fflush(stdout);
#endif
                        if ( candidateMate[v-StartIndex] == u ) { //Only if pointing to the matched vertex
                            //Start: PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
                            //Start: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
                            adj11 = verLocPtr[v-StartIndex];
                            adj12 = verLocPtr[v-StartIndex+1];
                            w = -1;
                            heaviestEdgeWt = MilanRealMin; //Assign the smallest Value possible first LDBL_MIN
                            for( k1 = adj11; k1 < adj12; k1++ ) {
                                if ( (verLocInd[k1]<StartIndex) || (verLocInd[k1]>EndIndex) ) { //Is it a ghost vertex?
                                    if(GMate[Ghost2LocalMap[verLocInd[k1]]] >= 0 )// Already matched
                                        continue;
                                }
                                else { //A local vertex
                                    if( Mate[verLocInd[k1]-StartIndex] >= 0 ) // Already matched
                                        continue;
                                }

                                if( (edgeLocWeight[k1] > heaviestEdgeWt) ||
                                    ((edgeLocWeight[k1] == heaviestEdgeWt)&&(w < verLocInd[k1])) ) {
                                    heaviestEdgeWt = edgeLocWeight[k1];
                                    w = verLocInd[k1];
                                }
                            } //End of for loop
                            candidateMate[v-StartIndex] = w;
                            //End: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
#ifdef PRINT_DEBUG_INFO_
                            cout<<"\n("<<myRank<<")"<<v<<" Points to: "<<w; fflush(stdout);
#endif
                            //If found a dominating edge:
                            if ( w >= 0 ) {
                                if ( (w < StartIndex) || (w > EndIndex) ) { //w is a ghost
                                    //Build the Message Packet:
                                    Message[0] = v; //LOCAL
                                    Message[1] = w; //GHOST
                                    Message[2] = REQUEST;  //TYPE
                                    //Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
                                    cout<<"\n("<<myRank<<")Sending a request message:";
                                    cout<<"\n("<<myRank<<")Ghost is "<<w<<" Owner is: "<<findOwnerOfGhost(w, verDistance, myRank, numProcs);
                                    fflush(stdout);
#endif
                                    ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs); assert(ghostOwner != -1); assert(ghostOwner != myRank);
                                    MPI_Bsend(&Message[0], 3, TypeMap<MilanLongInt>(), ghostOwner, ComputeTag, comm);
                                    msgInd++; msgActual++;
                                    if ( candidateMate[NLVer+Ghost2LocalMap[w]] == v ) {
                                        Mate[v-StartIndex] = w; //v is local
                                        GMate[Ghost2LocalMap[w]] = v; //w is ghost
                                        //Q.push_back(u);
                                        U.push_back(v);
                                        U.push_back(w);
                                        myCard++;
#ifdef PRINT_DEBUG_INFO_
                                        cout<<"\n("<<myRank<<")MATCH: ("<<v<<","<<w<<") "; fflush(stdout);
#endif
                                        //Decrement the counter:
                                        //Start: PARALLEL_PROCESS_CROSS_EDGE_B(v,w)
                                        if ( Counter[Ghost2LocalMap[w]] > 0 ) {
                                            Counter[Ghost2LocalMap[w]] = Counter[Ghost2LocalMap[w]] - 1; //Decrement
                                            if ( Counter[Ghost2LocalMap[w]] == 0 ) {
                                                S--; //Decrement S
#ifdef PRINT_DEBUG_INFO_
                                                cout<<"\n("<<myRank<<")Decrementing S: Ghost vertex "<<w<<" has received all its messages";
                                                fflush(stdout);
#endif
                                            }
                                        } //End of if Counter[w] > 0
                                        //End: PARALLEL_PROCESS_CROSS_EDGE_B(v,w)
                                    } //End of if CandidateMate[w] = v
                                } //End of if a Ghost Vertex
                                else { //w is a local vertex
                                    if ( candidateMate[w-StartIndex] == v )  {
                                        Mate[v-StartIndex] = w; //v is local
                                        Mate[w-StartIndex] = v; //w is local
                                        //Q.push_back(u);
                                        U.push_back(v);
                                        U.push_back(w);
                                        myCard++;
#ifdef PRINT_DEBUG_INFO_
                                        cout<<"\n("<<myRank<<")MATCH: ("<<v<<","<<w<<") "; fflush(stdout);
#endif
                                    } //End of if(CandidateMate(w) = v
                                } //End of Else
                            } //End of if(w >=0)
                            else { //no dominating edge found: w == -1
                                adj11 = verLocPtr[v-StartIndex];
                                adj12 = verLocPtr[v-StartIndex+1];
                                for( k1 = adj11; k1 < adj12; k1++ ) {
                                    w = verLocInd[k1];
                                    if ( (w < StartIndex) || (w > EndIndex) ) { //A ghost
                                        //Build the Message Packet:
                                        Message[0] = v;			 //LOCAL
                                        Message[1] = w;            //GHOST
                                        Message[2] = FAILURE;      //TYPE
                                        //Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
                                        cout<<"\n("<<myRank<<")Sending a failure message: ";
                                        cout<<"\n("<<myRank<<")Ghost is "<<w<<" Owner is: "<<findOwnerOfGhost(w, verDistance, myRank, numProcs);
                                        fflush(stdout);
#endif
                                        ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs); assert(ghostOwner != -1); assert(ghostOwner != myRank);
                                        MPI_Bsend(&Message[0], 3, TypeMap<MilanLongInt>(), ghostOwner, ComputeTag, comm);
                                        msgInd++; msgActual++;
                                    } //End of if(GHOST)
                                } //End of for loop
                            } // End of Else: w == -1
                            //End:   PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
                        } //End of If (candidateMate[v-StartIndex] == u)
                    } //End of if ( (v >= StartIndex) && (v <= EndIndex) ) //If Local Vertex:
                    else { //Neighbor v is a ghost vertex
                        if ( candidateMate[NLVer+Ghost2LocalMap[v]] == u )
                            candidateMate[NLVer+Ghost2LocalMap[v]] = -1;
                        if ( v != Mate[u-StartIndex] ) { //u is a local vertex
                            //Build the Message Packet:
                            Message[0] = u; //LOCAL
                            Message[1] = v; //GHOST
                            Message[2] = SUCCESS;  //TYPE
                            //Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
                            cout<<"\n("<<myRank<<")Sending a success message: ";
                            cout<<"\n("<<myRank<<")Ghost is "<<v<<" Owner is: "<<findOwnerOfGhost(v, verDistance, myRank, numProcs);
                            fflush(stdout);
#endif
                            ghostOwner = findOwnerOfGhost(v, verDistance, myRank, numProcs); assert(ghostOwner != -1); assert(ghostOwner != myRank);
                            MPI_Bsend(&Message[0], 3, TypeMap<MilanLongInt>(), ghostOwner, ComputeTag, comm);
                            msgInd++;
                            msgActual++;
#ifdef DEBUG_GHOST_
                            if ((u<StartIndex) || (u>EndIndex)) {
			      cout<<"\n("<<myRank<<") "<<__LINE__<<" From Send: should not happen: u= "<<u<<" v= "<<v<<
				" StartIndex "<<StartIndex<<" EndIndex "<<EndIndex<<endl;
			      fflush(stdout);
			    }
#endif

                        } //End of If( v != Mate[u] )
                    } //End of Else //A Ghost Vertex
                } //End of For Loop adj(u)
            } //End of if ( (u >= StartIndex) && (u <= EndIndex) ) //Process Only If a Local Vertex
        } //End of while ( /*!Q.empty()*/ !U.empty() )
        ///////////////////////// END OF PROCESS MATCHED VERTICES /////////////////////////

        //// BREAK IF NO MESSAGES EXPECTED /////////
#ifdef PRINT_DEBUG_INFO_
        cout<<"\n("<<myRank<<")Deciding whether to break: S= "<<S<<endl;
#endif

        if ( S == 0 ) {
#ifdef DEBUG_HANG_
            cout<<"\n("<<myRank<<") Breaking out" <<endl; fflush(stdout);
#endif
            break;
        }
        ///////////////////////////////////////////////////////////////////////////////////
        /////////////////////////// PROCESS MESSAGES //////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////
        /*
         RECEIVE message ( u, v, message_type );
         // u is a GHOST vertex ... v is a LOCAL vertex
         */
#ifdef PRINT_DEBUG_INFO_
        cout<<"\n("<<myRank<<"=========================************==============================="<<endl; fflush(stdout);
        fflush(stdout);
#endif
#ifdef PRINT_DEBUG_INFO_
        cout<<"\n("<<myRank<<")About to begin Message processing phase ... S="<<S<<endl;
        fflush(stdout);
#endif
#ifdef PRINT_DEBUG_INFO_
        cout<<"\n("<<myRank<<"=========================************==============================="<<endl; fflush(stdout);
        fflush(stdout);
#endif
        //BLOCKING RECEIVE:
#ifdef PRINT_DEBUG_INFO_
        cout<<"\n("<<myRank<<" Waiting for blocking receive..."<<endl; fflush(stdout);
        fflush(stdout);
#endif
        error_codeC = MPI_Recv(&Message[0], 3, TypeMap<MilanLongInt>(), MPI_ANY_SOURCE, ComputeTag, comm, &computeStatus);
        if (error_codeC != MPI_SUCCESS ) {
            MPI_Error_string(error_codeC, error_message, &message_length);
            cout<<"\n*Error in call to MPI_Receive on Slave: "<<error_message<<"\n"; fflush(stdout);
        }
        Sender = computeStatus.MPI_SOURCE;
#ifdef PRINT_DEBUG_INFO_
        cout<<"\n("<<myRank<<")Received message from Process "<<Sender<<" Type= "<<Message[2]<<endl;
        fflush(stdout);
#endif
        //If the Message Type is a size indicator, then receive the bigger message.
        if ( Message[2] == SIZEINFO ) {
#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")Received bundled message from Process "<<Sender<<" Size= "<<Message[0]<<endl;
            fflush(stdout);
#endif
            bundleSize = Message[0]; //#of integers in the message
            //Build the Message Buffer:
            if (!ReceiveBuffer.empty())
                ReceiveBuffer.clear(); //Empty it out first
            ReceiveBuffer.resize(bundleSize, -1); //Initialize
#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")Message Bundle Before: "<<endl;
            for (i=0; i<bundleSize; i++)
                cout<<ReceiveBuffer[i]<<",";
            cout<<endl;
            fflush(stdout);
#endif
            //Receive the message
            error_codeC = MPI_Recv(&ReceiveBuffer[0], bundleSize, TypeMap<MilanLongInt>(), Sender, BundleTag, comm, &computeStatus);
            if (error_codeC != MPI_SUCCESS ) {
                MPI_Error_string(error_codeC, error_message, &message_length);
                cout<<"\n*Error in call to MPI_Receive on processor "<<myRank<<" Error: "<<error_message<<"\n"; fflush(stdout);
            }
#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")Message Bundle After: "<<endl;
            for (i=0; i<bundleSize; i++)
                cout<<ReceiveBuffer[i]<<",";
            cout<<endl;
            fflush(stdout);
#endif
        }
        else { //Just a single message:
#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")Received regular message from Process "<<Sender<<" u= "<<Message[0]<<" v= "<<Message[1]<<endl;
            fflush(stdout);
#endif
            //Add the current message to Queue:
            bundleSize = 3; //#of integers in the message
            //Build the Message Buffer:
            if (!ReceiveBuffer.empty())
                ReceiveBuffer.clear(); //Empty it out first
            ReceiveBuffer.resize(bundleSize, -1); //Initialize

            ReceiveBuffer[0]=Message[0]; //u
            ReceiveBuffer[1]=Message[1]; //v
            ReceiveBuffer[2]=Message[2]; //message_type
        }
        bundleCounter = 0;
        while ( bundleCounter < bundleSize ) {
            u = ReceiveBuffer[bundleCounter]; //GHOST
            bundleCounter++;
            v = ReceiveBuffer[bundleCounter]; //LOCAL
            bundleCounter++;
            message_type = ReceiveBuffer[bundleCounter]; //TYPE
            bundleCounter++;
#ifdef DEBUG_GHOST_
            if ((v<StartIndex) || (v>EndIndex)) {
	      cout<<"\n("<<myRank<<") From ReceiveBuffer: This should not happen: u= "<<u<<" v= "<<v<<" Type= "<<message_type<<
		" StartIndex "<<StartIndex<<" EndIndex "<<EndIndex<<endl;
	      fflush(stdout);
	    }
#endif
#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")Processing message: u= "<<u<<" v= "<<v<<" Type= "<<message_type<<endl;
            fflush(stdout);
#endif
            // CASE I: REQUEST
            if ( message_type == REQUEST ) {
#ifdef PRINT_DEBUG_INFO_
                cout<<"\n("<<myRank<<")Message type is REQUEST"<<endl; fflush(stdout);
#endif
#ifdef DEBUG_GHOST_
                if ((v<0)||(v<StartIndex) || ((v-StartIndex)>NLVer)) {
		  cout<<"\n("<<myRank<<") case 1 Bad address "<<v<<" "<<StartIndex<<" "<<v-StartIndex<<" "<<NLVer<<endl; fflush(stdout);
		}

#endif
                if ( Mate[v-StartIndex] == -1 ) { //Process only if not already matched  (v is local)
                    candidateMate[NLVer+Ghost2LocalMap[u]] = v;  //Set CandidateMate for the ghost
                    if ( candidateMate[v-StartIndex] == u ) {
                        GMate[Ghost2LocalMap[u]] = v; //u is ghost
                        Mate[v-StartIndex] = u; //v is local
                        //Q.push_back(u);
                        U.push_back(v);
                        U.push_back(u);
                        myCard++;
#ifdef PRINT_DEBUG_INFO_
                        cout<<"\n("<<myRank<<")MATCH: ("<<v<<","<<u<<") "<<endl; fflush(stdout);
#endif
                        //Start: PARALLEL_PROCESS_CROSS_EDGE_B(v,u)
                        if ( Counter[Ghost2LocalMap[u]] > 0 ) {
                            Counter[Ghost2LocalMap[u]] = Counter[Ghost2LocalMap[u]] - 1; //Decrement
                            if ( Counter[Ghost2LocalMap[u]] == 0 ) {
                                S--; //Decrement S
#ifdef PRINT_DEBUG_INFO_
                                cout<<"\n("<<myRank<<")Decrementing S: Ghost vertex "<<u<<" has received all its messages"<<endl;
                                fflush(stdout);
#endif
                            }
                        } //End of if Counter[w] > 0
                        //End: PARALLEL_PROCESS_CROSS_EDGE_B(v,u)
                    } //End of if ( candidateMate[v-StartIndex] == u )e
                } //End of if ( Mate[v] == -1 )
            } //End of REQUEST
            else {   //CASE II: SUCCESS
                if ( message_type == SUCCESS ) {
#ifdef PRINT_DEBUG_INFO_
                    cout<<"\n("<<myRank<<")Message type is SUCCESS"<<endl; fflush(stdout);
#endif
                    //Start: PARALLEL_PROCESS_CROSS_EDGE_B(v,u)
                    GMate[Ghost2LocalMap[u]] = EndIndex+1; //Set a Dummy Mate to make sure that we do not (u is a ghost)
                    //process it again
                    if ( Counter[Ghost2LocalMap[u]] > 0 ) {
                        Counter[Ghost2LocalMap[u]] = Counter[Ghost2LocalMap[u]] - 1; //Decrement
                        if ( Counter[Ghost2LocalMap[u]] == 0 ) {
                            S--; //Decrement S
#ifdef PRINT_DEBUG_INFO_
                            cout<<"\n("<<myRank<<")Decrementing S: Ghost vertex "<<u<<" has received all its messages";
                            fflush(stdout);
#endif
                        }
                    } //End of if Counter[w] > 0
                    //End: PARALLEL_PROCESS_CROSS_EDGE_B(v,u)
#ifdef DEBUG_GHOST_
                    if ((v<0)||(v<StartIndex) || ((v-StartIndex)>NLVer)) {
		      cout<<"\n("<<myRank<<") case 2  Bad address "<<v<<" "<<StartIndex<<" "<<v-StartIndex<<" "<<NLVer<<endl; fflush(stdout);
		    }
#endif
                    if ( Mate[v-StartIndex] == -1 ) { //Process only if not already matched ( v is local)
                        if ( candidateMate[v-StartIndex] == u ) {
                            //Start: PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
                            //Start: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
                            adj11 = verLocPtr[v-StartIndex];
                            adj12 = verLocPtr[v-StartIndex+1];
                            w = -1;
                            heaviestEdgeWt = MilanRealMin; //Assign the smallest Value possible first LDBL_MIN
                            for( k1 = adj11; k1 < adj12; k1++ ) {
                                if ( (verLocInd[k1]<StartIndex) || (verLocInd[k1]>EndIndex) ) { //Is it a ghost vertex?
                                    if(GMate[Ghost2LocalMap[verLocInd[k1]]] >= 0 )// Already matched
                                        continue;
                                }
                                else { //A local vertex
                                    if( Mate[verLocInd[k1]-StartIndex] >= 0 ) // Already matched
                                        continue;
                                }

                                if( (edgeLocWeight[k1] > heaviestEdgeWt) ||
                                    ((edgeLocWeight[k1] == heaviestEdgeWt)&&(w < verLocInd[k1])) ) {
                                    heaviestEdgeWt = edgeLocWeight[k1];
                                    w = verLocInd[k1];
                                }
                            } //End of for loop
                            candidateMate[v-StartIndex] = w;
                            //End: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
#ifdef PRINT_DEBUG_INFO_
                            cout<<"\n("<<myRank<<")"<<v<<" Points to: "<<w<<endl; fflush(stdout);
#endif
                            //If found a dominating edge:
                            if ( w >= 0 ) {
                                if ( (w < StartIndex) || (w > EndIndex) ) { //w is a ghost
                                    //Build the Message Packet:
                                    Message[0] = v; //LOCAL
                                    Message[1] = w; //GHOST
                                    Message[2] = REQUEST;  //TYPE
                                    //Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
                                    cout<<"\n("<<myRank<<")Sending a request message: ";
                                    cout<<"\n("<<myRank<<")Ghost is "<<w<<" Owner is: "<<findOwnerOfGhost(w, verDistance, myRank, numProcs)<<endl;
                                    fflush(stdout);
#endif
                                    ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs); assert(ghostOwner != -1); assert(ghostOwner != myRank);
                                    MPI_Bsend(&Message[0], 3, TypeMap<MilanLongInt>(), ghostOwner, ComputeTag, comm);
                                    msgInd++; msgActual++;
                                    if ( candidateMate[NLVer+Ghost2LocalMap[w]] == v ) {
                                        Mate[v-StartIndex] = w; //v is local
                                        GMate[Ghost2LocalMap[w]] = v; //w is ghost
                                        //Q.push_back(u);
                                        U.push_back(v);
                                        U.push_back(w);
                                        myCard++;
#ifdef PRINT_DEBUG_INFO_
                                        cout<<"\n("<<myRank<<")MATCH: ("<<v<<","<<w<<") "<<endl; fflush(stdout);
#endif
                                        //Decrement the counter:
                                        //Start: PARALLEL_PROCESS_CROSS_EDGE_B(v,w)
                                        if ( Counter[Ghost2LocalMap[w]] > 0 ) {
                                            Counter[Ghost2LocalMap[w]] = Counter[Ghost2LocalMap[w]] - 1; //Decrement
                                            if ( Counter[Ghost2LocalMap[w]] == 0 ) {
                                                S--; //Decrement S
#ifdef PRINT_DEBUG_INFO_
                                                cout<<"\n("<<myRank<<")Decrementing S: Ghost vertex "<<w<<" has received all its messages";
                                                fflush(stdout);
#endif
                                            }
                                        } //End of if Counter[w] > 0
                                        //End: PARALLEL_PROCESS_CROSS_EDGE_B(v,w)
                                    } //End of if CandidateMate[w] = v
                                } //End of if a Ghost Vertex
                                else { //w is a local vertex
                                    if ( candidateMate[w-StartIndex] == v ) {
                                        Mate[v-StartIndex] = w; //v is local
                                        Mate[w-StartIndex] = v; //w is local
                                        //Q.push_back(u);
                                        U.push_back(v);
                                        U.push_back(w);
                                        myCard++;
#ifdef PRINT_DEBUG_INFO_
                                        cout<<"\n("<<myRank<<")MATCH: ("<<v<<","<<w<<") "<<endl; fflush(stdout);
#endif
                                    } //End of if(CandidateMate(w) = v
                                } //End of Else
                            } //End of if(w >=0)
                            else { //No dominant edge found
                                adj11 = verLocPtr[v-StartIndex];
                                adj12 = verLocPtr[v-StartIndex+1];
                                for( k1 = adj11; k1 < adj12; k1++ ) {
                                    w = verLocInd[k1];
                                    if ( (w < StartIndex) || (w > EndIndex) ) { //A ghost
                                        //Build the Message Packet:
                                        Message[0] = v;			 //LOCAL
                                        Message[1] = w;            //GHOST
                                        Message[2] = FAILURE;      //TYPE
                                        //Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
                                        cout<<"\n("<<myRank<<")Sending a failure message: ";
                                        cout<<"\n("<<myRank<<")Ghost is "<<w<<" Owner is: "<<findOwnerOfGhost(w, verDistance, myRank, numProcs)<<endl;
                                        fflush(stdout);
#endif
                                        //MPI_Bsend(&Message[0], 3, MilanMpiLongInt, findOwnerOfGhost(w, verDistance, myRank, numProcs),
                                        ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs); assert(ghostOwner != -1); assert(ghostOwner != myRank);
                                        MPI_Bsend(&Message[0], 3, TypeMap<MilanLongInt>(), ghostOwner, ComputeTag, comm);
                                        msgInd++; msgActual++;
                                    } //End of if(GHOST)
                                } //End of for loop
                            } // End of Else: w == -1
                            //End:   PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
                        } //End of if ( candidateMate[v-StartIndex] == u )
                    } //End of if ( Mate[v] == -1 )
                } //End of if ( message_type == SUCCESS )
                else { //CASE III: FAILURE
#ifdef PRINT_DEBUG_INFO_
                    cout<<"\n("<<myRank<<")Message type is FAILURE"<<endl; fflush(stdout);
#endif
                    //Start: PARALLEL_PROCESS_CROSS_EDGE_B(v,u)
                    GMate[Ghost2LocalMap[u]] = EndIndex+1; //Set a Dummy Mate to make sure that we do not (u is a ghost)
                    //process it again
                    if ( Counter[Ghost2LocalMap[u]] > 0 ) {
                        Counter[Ghost2LocalMap[u]] = Counter[Ghost2LocalMap[u]] - 1; //Decrement
                        if ( Counter[Ghost2LocalMap[u]] == 0 ) {
                            S--; //Decrement S
#ifdef PRINT_DEBUG_INFO_
                            cout<<"\n("<<myRank<<")Decrementing S: Ghost vertex "<<u<<" has received all its messages";
                            fflush(stdout);
#endif
                        }
                    } //End of if Counter[w] > 0
                    //End: PARALLEL_PROCESS_CROSS_EDGE_B(v,u)
                } //End of else: CASE III
            } //End of else: CASE I
        } //End of if (!MsgQ.empty())
        ///////////////////////// END OF PROCESS MESSAGES /////////////////////////////////
#ifdef PRINT_DEBUG_INFO_
        cout<<"\n("<<myRank<<")Finished Message processing phase: S= "<<S; fflush(stdout);
        cout<<"\n("<<myRank<<")** SENT     : ACTUAL= "<<msgActual; fflush(stdout);
        cout<<"\n("<<myRank<<")** SENT     : INDIVIDUAL= "<<msgInd<<endl; fflush(stdout);
#endif
    } //End of while (true)

#ifdef PRINT_DEBUG_INFO_
    cout<<"\n("<<myRank<<") Waitall= "<<endl; fflush(stdout);
#endif
#ifdef DEBUG_HANG_
    cout<<"\n("<<myRank<<") Waitall " <<endl; fflush(stdout);
#endif
    //MPI_Barrier(comm);
    //Cleanup Phase
    MPI_Waitall(MessageIndex, &SRequest[0], &SStatus[0]);
    //MPI_Buffer_attach(&Buffer, BufferSize); //Attach the Buffer
    if ( BufferSize > 0 ) {
        MPI_Buffer_detach(&Buffer, &BufferSize); //Detach the Buffer
        free(Buffer); //Free the memory that was allocated
    }
    finishTime = MPI_Wtime();
    *ph2_time = finishTime-startTime; //Time taken for Phase-2
    *ph2_card = myCard ; //Cardinality at the end of Phase-2

#ifdef PRINT_DEBUG_INFO_
    cout<<"\n("<<myRank<<")End of function to compute matching: "<<endl; fflush(stdout);
    cout<<"\n("<<myRank<<")myCardinality: "<<myCard<<endl; fflush(stdout);
    cout<<"\n("<<myRank<<")Matching took "<<finishTime-startTime<<"seconds"<<endl; fflush(stdout);
    cout<<"\n("<<myRank<<")** Getting out of the matching function **"<<endl; fflush(stdout);
#endif
#ifdef PRINT_DEBUG_INFO_
    cout<<"\n("<<myRank<<") Number of Ghost edges = "<<numGhostEdges;
    cout<<"\n("<<myRank<<") Total number of potential message X 2 = "<<numGhostEdges*2;
    cout<<"\n("<<myRank<<") Number messages bundled = "<<NumMessagesBundled;
    cout<<"\n("<<myRank<<") Total Individual Messages sent = "<< msgInd;
    if (msgInd>0) {
      cout<<"\n("<<myRank<<") Percentage of messages bundled = "<<((double)NumMessagesBundled/(double)(msgInd))*100.0<<"% \n";
    }
    fflush(stdout);
#endif

    *msgActualSent = msgActual;
    *msgIndSent = msgInd;
    if (msgInd > 0) {
        *msgPercent = ((double)NumMessagesBundled/(double)(msgInd))*100.0;
    } else {
        *msgPercent = 0;
    }

#ifdef DEBUG_HANG_
    if (myRank == 0) cout<<"\n("<<myRank<<") Done" <<endl; fflush(stdout);
#endif
    //MPI_Barrier(comm);
}
//End of algoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMate
#endif

#endif