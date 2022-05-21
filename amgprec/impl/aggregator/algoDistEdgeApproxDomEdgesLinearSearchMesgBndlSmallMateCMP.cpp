#include "MatchBoxPC.h"
#include <omp.h>
#include <stdio.h>
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
    //Get the iterators for the graph:
    //vector<MilanLongInt>::iterator verLocPtr  = inputSubGraph.getVerPtr_b();
    //vector<MilanLongInt>::iterator verLocInd  = inputSubGraph.getVerInd_b();
    //vector<MilanReal>::iterator edgeLocWeight = inputSubGraph.getEdgeWt_b();

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

    MilanLongInt NumMessagesBundled;
    MilanInt ghostOwner; // Changed by Fabio to be an integer, addresses needs to be integers!
    //vector<MilanLongInt> candidateMate;
    MilanLongInt* candidateMate = new MilanLongInt[1];
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
    // index that starts with zero to |Vg|  - 1
    map<MilanLongInt, MilanLongInt>::iterator storedAlready;
    vector <MilanLongInt> Counter;  //Store the edge count for each ghost vertex
    MilanLongInt numGhostVertices = 0, numGhostEdges = 0, insertMe = 0; //Number of Ghost vertices

#ifdef PRINT_DEBUG_INFO_
    cout<<"\n("<<myRank<<")About to compute Ghost Vertices..."; fflush(stdout);
#endif
#ifdef DEBUG_HANG_
    if (myRank == 0)     cout<<"\n("<<myRank<<")About to compute Ghost Vertices..."; fflush(stdout);
#endif

    /*
     * OMP Ghost2LocalInitialization
     * The cycle analyzes all the edges and when finds a ghost edge
     * puts it in the Ghost2LocalMap.
     * A critical region is needed when inserting data in the map.
     *
     * Despite the critical region it is still productive to
     * parallelize this for because the critical region is exeuted
     * only when a ghost edge is found and ghost edges are a minority.
     */

    //Define Adjacency Lists for Ghost Vertices:
    //cout<<"Building Ghost data structures ... \n\n";
    vector <MilanLongInt> verGhostPtr, verGhostInd, tempCounter;
    //Mate array for ghost vertices:
    vector <MilanLongInt> GMate;  //Proportional to the number of ghost vertices
    MilanLongInt S;
    staticQueue U;
#ifdef TIME_TRACKER
    double Ghost2LocalInitialization = MPI_Wtime();
#endif

#pragma omp parallel private(insertMe, k, k1, adj1, adj2, adj11, adj12, heaviestEdgeWt, w, ghostOwner) firstprivate(StartIndex, EndIndex) default(shared) num_threads(4)
    {

        // TODO comments about the reduction

#pragma omp for reduction(+ : numGhostEdges)
        for (i = 0; i < NLEdge; i++) { //O(m) - Each edge stored twice
            insertMe = verLocInd[i];
            //cout<<"InsertMe on Process "<<myRank<<" is: "<<insertMe<<endl;
            if ((insertMe < StartIndex) || (insertMe > EndIndex)) { //Find a ghost
                numGhostEdges++;
#pragma omp critical
                {
                    storedAlready = Ghost2LocalMap.find(insertMe);
                    if (storedAlready != Ghost2LocalMap.end()) { //Has already been added
                        //cout<<"Process "<<myRank<<" found: "<<storedAlready->first<<" - "<<storedAlready->second<<endl;
                        Counter[storedAlready->second]++; //Increment the counter
                    } else { //Insert an entry for the ghost:
                        //cout<<"Process "<<myRank<<" * New insert:  Key="<<insertMe<< " : Value="<<numGhostVertices<<endl;
                        Ghost2LocalMap[insertMe] = numGhostVertices; //Add a map entry
                        Counter.push_back(1); //Initialize the counter
                        numGhostVertices++;  //Increment the number of ghost vertices
                    } //End of else()
                }
            } //End of if ( (insertMe < StartIndex) || (insertMe > EndIndex) )
        } //End of for(ghost vertices)

#pragma omp single
        {
            //numGhostEdges = atomicNumGhostEdges;
#ifdef TIME_TRACKER
            Ghost2LocalInitialization = MPI_Wtime() - Ghost2LocalInitialization;
            fprintf(stderr, "Ghost2LocalInitialization time: %f\n", Ghost2LocalInitialization);
#endif

#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")NGhosts:" << numGhostVertices << " GhostEdges: "<<numGhostEdges;
            if (!Ghost2LocalMap.empty()) {
                cout<<"\n("<<myRank<<")Final Map : on process ";
                cout<<"\n("<<myRank<<")Key \t Value \t Counter \n"; fflush(stdout);
                storedAlready = Ghost2LocalMap.begin();
                do {
                    cout<<storedAlready->second<<" - "<<storedAlready->first<<" : "<<Counter[storedAlready->second]<<endl;
                    fflush(stdout);
                    storedAlready++;
                } while ( storedAlready != Ghost2LocalMap.end() );
            }
#endif

            //Initialize adjacency Lists for Ghost Vertices:
            try {
                verGhostPtr.reserve(numGhostVertices + 1); //Pointer Vector
                tempCounter.reserve(numGhostVertices); //Pointer Vector
                verGhostInd.reserve(numGhostEdges); //Index Vector
                GMate.reserve(numGhostVertices); //Ghost Mate Vector
            } catch (length_error) {
                cout << "Error in function algoDistEdgeApproxDominatingEdgesLinearSearch: \n";
                cout << "Not enough memory to allocate the internal variables \n";
                exit(1);
            }
            //Initialize the Vectors:
            verGhostPtr.resize(numGhostVertices + 1, 0); //Pointer Vector
            tempCounter.resize(numGhostVertices, 0); //Temporary Counter
            verGhostInd.resize(numGhostEdges, -1); //Index Vector
            GMate.resize(numGhostVertices, -1); //Temporary Counter
            verGhostPtr[0] = 0; //The first value
#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")Ghost Vertex Pointer: "; fflush(stdout);
#endif

#ifdef TIME_TRACKER
            double verGhostPtrInitialization = MPI_Wtime();
#endif

        }
        /*
         * OMP verGhostPtrInitialization
         *
         */
#pragma omp for nowait
        for (i = 0; i < numGhostVertices; i++) { //O(|Ghost Vertices|)
            verGhostPtr[i + 1] = verGhostPtr[i] + Counter[i];
#ifdef PRINT_DEBUG_INFO_
            cout<<verGhostPtr[i]<<"\t"; fflush(stdout);
#endif
        }

#ifdef TIME_TRACKER
        verGhostPtrInitialization = MPI_Wtime() - verGhostPtrInitialization;
        fprintf(stderr, "verGhostPtrInitialization time: %f\n", verGhostPtrInitialization);
#endif

#ifdef PRINT_DEBUG_INFO_
        if ( numGhostVertices > 0 )
            cout<<verGhostPtr[numGhostVertices]<<"\n";
        fflush(stdout);
#endif

        /*
         * OMP verGhostIndInitialization
         *
         * In this cycle the verGhostInd is initialized
         * with the datas related to ghost edges.
         * The check to see if a node is a ghost node is
         * executed in paralle and when a ghost node
         * is found a critical region is started.
         *
         * Despite the critical region it's still useful to
         * parallelize the for cause the ghost nodes
         * are a minority hence the critical region is executed
         * few times.
         */

#ifdef TIME_TRACKER
        double verGhostIndInitialization = MPI_Wtime();
#endif

#pragma omp for nowait
        for (v = 0; v < NLVer; v++) {
            adj1 = verLocPtr[v];   //Vertex Pointer
            adj2 = verLocPtr[v + 1];
            for (k = adj1; k < adj2; k++) {
                w = verLocInd[k]; //Get the adjacent vertex
                if ((w < StartIndex) || (w > EndIndex)) { //Find a ghost
#pragma omp critical
                    {
                        insertMe = verGhostPtr[Ghost2LocalMap[w]] + tempCounter[Ghost2LocalMap[w]]; //Where to insert
                        verGhostInd[insertMe] = v + StartIndex; //Add the adjacency
                        tempCounter[Ghost2LocalMap[w]]++; //Increment the counter
                    }
                } //End of if((w < StartIndex) || (w > EndIndex))
            } //End of for(k)
        } //End of for (v)

#pragma omp single
        {

#ifdef TIME_TRACKER
            verGhostIndInitialization = MPI_Wtime() - verGhostIndInitialization;
            fprintf(stderr, "verGhostIndInitialization time: %f\n", verGhostIndInitialization);
#endif

#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")Ghost Vertex Index: ";
            for ( v=0; v < numGhostEdges; v++ )
                cout<<verGhostInd[v]<<"\t";
            cout<<endl; fflush(stdout);
#endif


            Message.resize(3, -1);
            message_type = 0;
            NumMessagesBundled = 0;
            ghostOwner = 0;
            try {
                QLocalVtx.reserve(numGhostEdges); //Local Vertex
                QGhostVtx.reserve(numGhostEdges); //Ghost Vertex
                QMsgType.reserve(numGhostEdges); //Message Type (Request/Failure)
                QOwner.reserve(numGhostEdges); //Owner of the ghost: COmpute once and use later
            } catch (length_error) {
                cout << "Error in function algoDistEdgeApproxDominatingEdgesMessageBundling: \n";
                cout << "Not enough memory to allocate the internal variables \n";
                exit(1);
            }

#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<")Allocating CandidateMate.. "; fflush(stdout);
#endif
            //Allocate Data Structures:
            /*
             * candidateMate was a vector and has been replaced with a raw array
             * there is no point in using the vector (or maybe there is???)
             * so I replaced it with an array wich is slightly faster
             */
            delete[] candidateMate;
            candidateMate = new MilanLongInt[NLVer + numGhostVertices];

            /*
             * Create the Queue Data Structure for the Dominating Set
             *
             * I had to declare the staticuQueue U before the parallel region
             * to have it in the correct scope. Since we can't chane the dimension
             * of a staticQueue I had to destroy the previous object and instantiate
             * a new one of the correct size.
             */
            U.~staticQueue();
            new(&U) staticQueue(NLVer + numGhostVertices);

#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<"=========================************==============================="<<endl; fflush(stdout);
            fflush(stdout);
#endif
            //MPI_Barrier(comm);
            finishTime = MPI_Wtime();
            *ph0_time = finishTime - startTime; //Time taken for Phase-0: Initialization
#ifdef PRINT_DEBUG_INFO_
            cout<<"\n("<<myRank<<") Setup Time :"<< *ph0_time <<endl; fflush(stdout);
            fflush(stdout);
#endif
#ifdef DEBUG_HANG_
            if (myRank == 0) cout<<"\n("<<myRank<<") Setup Time :"<< *ph0_time <<endl; fflush(stdout);
#endif
            startTime = MPI_Wtime();
            /////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////// INITIALIZATION /////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////
            //Compute the Initial Matching Set:

            S = numGhostVertices; //Initialize S with number of Ghost Vertices

            /*
             * OMP PARALLEL_COMPUTE_CANDIDATE_MATE_B
             * It is actually not possible to parallelize this cycle
             * as it is.
             *
             * TODO think how it could be parallelizable
             */

            for ( v=0; v < NLVer; v++ ) {
#ifdef PRINT_DEBUG_INFO_
                cout<<"\n("<<myRank<<")Processing: "<<v+StartIndex<<endl; fflush(stdout);
#endif
                //Start: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)
                adj1 = verLocPtr[v];
                adj2 = verLocPtr[v + 1];
                w = -1;
                heaviestEdgeWt = MilanRealMin; //Assign the smallest Value possible first LDBL_MIN
                for (k = adj1; k < adj2; k++) {
                    if (isAlreadyMatched(k, verLocInd, StartIndex, EndIndex, GMate, Mate, Ghost2LocalMap)) continue;

                    if ((edgeLocWeight[k] > heaviestEdgeWt) ||
                        ((edgeLocWeight[k] == heaviestEdgeWt) && (w < verLocInd[k]))) {
                        heaviestEdgeWt = edgeLocWeight[k];
                        w = verLocInd[k];
                    }
                } //End of for loop
                //printf("Compare %ld, %ld\n", w, firstComputeCandidateMate(verLocPtr[v], verLocPtr[v + 1], verLocInd, edgeLocWeight));
                candidateMate[v] = w;
                //End: PARALLEL_COMPUTE_CANDIDATE_MATE_B(v)

                //Start: PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)

#ifdef PRINT_DEBUG_INFO_
                cout<<"\n("<<myRank<<")Processing: "<<v+StartIndex<<endl; fflush(stdout);
#endif

#ifdef PRINT_DEBUG_INFO_
                cout<<"\n("<<myRank<<")"<<v+StartIndex<<" Points to: "<<w; fflush(stdout);
#endif
                //If found a dominating edge:
                if (w >= 0) {
                    myCard++;
                    if ((w < StartIndex) || (w > EndIndex)) { //w is a ghost vertex
                        //Build the Message Packet:
                        //Message[0] = v+StartIndex; //LOCAL
                        //Message[1] = w;            //GHOST
                        //Message[2] = REQUEST;      //TYPE
                        //Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
                        cout<<"\n("<<myRank<<")Sending a request message (291):";
                    cout<<"\n("<<myRank<<")Local is: "<<v+StartIndex<<" Ghost is "<<w<<" Owner is: "<< findOwnerOfGhost(w, verDistance, myRank, numProcs) <<endl;
                    fflush(stdout);
#endif
                        /* MPI_Bsend(&Message[0], 3, MPI_INT, inputSubGraph.findOwner(w),
                         ComputeTag, comm);*/
                        msgInd++;
                        NumMessagesBundled++;
                        ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
                        PCounter[ghostOwner]++;

                        QLocalVtx.push_back(v + StartIndex);
                        QGhostVtx.push_back(w);
                        QMsgType.push_back(REQUEST);
                        //ghostOwner = inputSubGraph.findOwner(w);
                        assert(ghostOwner != -1);
                        assert(ghostOwner != myRank);
                        QOwner.push_back(ghostOwner);

                        if (candidateMate[NLVer + Ghost2LocalMap[w]] == v + StartIndex) {

                            Mate[v] = w;
                            GMate[Ghost2LocalMap[w]] = v + StartIndex; //w is a Ghost
                            U.push_back(v + StartIndex);
                            U.push_back(w);

#ifdef PRINT_DEBUG_INFO_
                            cout<<"\n("<<myRank<<")MATCH: ("<<v+StartIndex<<","<<w<<")"; fflush(stdout);
#endif
                            //Decrement the counter:
                            //Start: PARALLEL_PROCESS_CROSS_EDGE_B(v)
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
                            //End: PARALLEL_PROCESS_CROSS_EDGE_B(v)
                        } //End of if CandidateMate[w] = v
                    } //End of if a Ghost Vertex
                    else { // w is a local vertex

                        if (candidateMate[w - StartIndex] == (v + StartIndex)) {

                            Mate[v] = w;  //v is local
                            Mate[w - StartIndex] = v + StartIndex; //w is local
                            //Q.push_back(u);
                            U.push_back(v + StartIndex);
                            U.push_back(w);

#ifdef PRINT_DEBUG_INFO_
                            cout<<"\n("<<myRank<<")MATCH: ("<<v+StartIndex<<","<<w<<") "; fflush(stdout);
#endif

                        } //End of if ( candidateMate[w-StartIndex] == (v+StartIndex) )
                    } //End of Else
                } //End of if(w >=0)
                else {
                    adj11 = verLocPtr[v];
                    adj12 = verLocPtr[v + 1];
                    for (k1 = adj11; k1 < adj12; k1++) {
                        w = verLocInd[k1];
                        if ((w < StartIndex) || (w > EndIndex)) { //A ghost
                            //Build the Message Packet:
                            //Message[0] = v+StartIndex; //LOCAL
                            //Message[1] = w;            //GHOST
                            //Message[2] = FAILURE;      //TYPE
                            //Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
                            cout<<"\n("<<myRank<<")Sending a failure message: ";
                        cout<<"\n("<<myRank<<")Ghost is "<<w<<" Owner is: "<<findOwnerOfGhost(w, verDistance, myRank, numProcs);
                        fflush(stdout);
#endif
                            /* MPI_Bsend(&Message[0], 3, MPI_INT, inputSubGraph.findOwner(w),
                             ComputeTag, comm); */
                            NumMessagesBundled++;
                            msgInd++;
                            ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
                            PCounter[ghostOwner]++;
#pragma omp critical
                            {
                                QLocalVtx.push_back(v + StartIndex);
                                QGhostVtx.push_back(w);
                                QMsgType.push_back(FAILURE);
                                //ghostOwner = inputSubGraph.findOwner(w);
                                assert(ghostOwner != -1);
                                assert(ghostOwner != myRank);
                                QOwner.push_back(ghostOwner);
                            }

                        } //End of if(GHOST)
                    } //End of for loop
                } // End of Else: w == -1
                //End:   PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
            } //End of for ( v=0; v < NLVer; v++ )

        } // end of single region
    } // end of parallel region

    tempCounter.clear(); //Do not need this any more


#ifdef PRINT_DEBUG_INFO_
    cout<<"\n("<<myRank<<"=========================************==============================="<<endl; fflush(stdout);
    fflush(stdout);
#endif
    ///////////////////////////////////////////////////////////////////////////////////
    /////////////////////////// PROCESS MATCHED VERTICES //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    while ( !U.empty() ) {
        u = U.pop_front(); //Get an element from the queue
#ifdef PRINT_DEBUG_INFO_
        cout<<"\n("<<myRank<<")u: "<<u; fflush(stdout);
#endif
        if ( (u >= StartIndex) && (u <= EndIndex) ) { //Process Only the Local Vertices
            //Get the Adjacency list for u
            adj1 = verLocPtr[u-StartIndex];  //Pointer
            adj2 = verLocPtr[u-StartIndex+1];
            for( k = adj1; k < adj2; k++ ) {
                v = verLocInd[k];
                if ( (v >= StartIndex) && (v <= EndIndex) ) { //If Local Vertex:
                    if ( (v<StartIndex) || (v>EndIndex) ) { //Is it a ghost vertex?
                        if(GMate[Ghost2LocalMap[v]] >= 0 )// Already matched
                            continue;
                    } else { //A local vertex
                        if( Mate[v-StartIndex] >= 0 ) // Already matched
                            continue;
                    } //End of else

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
                            } else { //A local vertex
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
                            if ( (w < StartIndex) || (w > EndIndex) ) { //A ghost
                                //Build the Message Packet:
                                //Message[0] = v; //LOCAL
                                //Message[1] = w; //GHOST
                                //Message[2] = REQUEST;  //TYPE
                                //Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
                                cout<<"\n("<<myRank<<")Sending a request message:";
                                cout<<"\n("<<myRank<<")Ghost is "<<w<<" Owner is: "<<findOwnerOfGhost(w, verDistance, myRank, numProcs);
#endif
                                /*MPI_Bsend(&Message[0], 3, MPI_INT, inputSubGraph.findOwner(w),
                                 ComputeTag, comm);*/
                                QLocalVtx.push_back(v);
                                QGhostVtx.push_back(w);
                                QMsgType.push_back(REQUEST);
                                //ghostOwner = inputSubGraph.findOwner(w);
                                ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs); assert(ghostOwner != -1); assert(ghostOwner != myRank);
                                QOwner.push_back(ghostOwner);
                                PCounter[ghostOwner]++;
                                NumMessagesBundled++;
                                msgInd++;
                                if ( candidateMate[NLVer+Ghost2LocalMap[w]] == v ) {
                                    Mate[v-StartIndex] = w;  //v is a local vertex
                                    GMate[Ghost2LocalMap[w]] = v;  //w is a ghost vertex
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
                                    Mate[v-StartIndex] = w;  //v is a local vertex
                                    Mate[w-StartIndex] = v;  //w is a local vertex
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
                        else {
                            adj11 = verLocPtr[v-StartIndex];
                            adj12 = verLocPtr[v-StartIndex+1];
                            for( k1 = adj11; k1 < adj12; k1++ ) {
                                w = verLocInd[k1];
                                if ( (w < StartIndex) || (w > EndIndex) ) { //A ghost
                                    //Build the Message Packet:
                                    //Message[0] = v;	     //LOCAL
                                    //Message[1] = w;            //GHOST
                                    //Message[2] = FAILURE;      //TYPE
                                    //Send a Request (Asynchronous)
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
                                    ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs); assert(ghostOwner != -1); assert(ghostOwner != myRank);
                                    QOwner.push_back(ghostOwner);
                                    PCounter[ghostOwner]++;
                                    NumMessagesBundled++;
                                    msgInd++;
                                } //End of if(GHOST)
                            } //End of for loop
                        } // End of Else: w == -1
                        //End:   PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
                    } //End of If (candidateMate[v-StartIndex] == u)
                } //End of if ( (v >= StartIndex) && (v <= EndIndex) ) //If Local Vertex:
                else { //Neighbor is a ghost vertex
                    if ( candidateMate[NLVer+Ghost2LocalMap[v]] == u )
                        candidateMate[NLVer+Ghost2LocalMap[v]] = -1;
                    if ( v != Mate[u-StartIndex] ) { //u is local
                        //Build the Message Packet:
                        //Message[0] = u; //LOCAL
                        //Message[1] = v; //GHOST
                        //Message[2] = SUCCESS;  //TYPE
                        //Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
                        cout<<"\n("<<myRank<<")Sending a success message: ";
                        cout<<"\n("<<myRank<<")Ghost is "<<v<<" Owner is: "<<findOwnerOfGhost(v, verDistance, myRank, numProcs)<<"\n"; fflush(stdout);
#endif
                        /* MPI_Bsend(&Message[0], 3, MPI_INT, inputSubGraph.findOwner(v),
                         ComputeTag, comm); */
                        QLocalVtx.push_back(u);
                        QGhostVtx.push_back(v);
                        QMsgType.push_back(SUCCESS);
                        //ghostOwner = inputSubGraph.findOwner(v);
                        ghostOwner = findOwnerOfGhost(v, verDistance, myRank, numProcs); assert(ghostOwner != -1); assert(ghostOwner != myRank);
                        QOwner.push_back(ghostOwner);
                        PCounter[ghostOwner]++;
                        NumMessagesBundled++;
                        msgInd++;
                    } //End of If( v != Mate[u] )
                } //End of Else //A Ghost Vertex
            } //End of For Loop adj(u)
        } //End of if ( (u >= StartIndex) && (u <= EndIndex) ) //Process Only If a Local Vertex
    } //End of while ( /*!Q.empty()*/ !U.empty() )
    ///////////////////////// END OF PROCESS MATCHED VERTICES /////////////////////////
#ifdef DEBUG_HANG_
    if (myRank == 0) cout<<"\n("<<myRank<<") Send Bundles" <<endl; fflush(stdout);
#endif
    /////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////// SEND BUNDLED MESSAGES /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    //Data structures for Bundled Messages:
    vector<MilanLongInt> PCumulative, PMessageBundle, PSizeInfoMessages;
    MilanLongInt myIndex=0;
    try {
        PMessageBundle.reserve(NumMessagesBundled*3); //Three integers per message
        PCumulative.reserve(numProcs+1); //Similar to Row Pointer vector in CSR data structure
        PSizeInfoMessages.reserve(numProcs*3); //Buffer to hold the Size info message packets
    } catch ( length_error ) {
        cout<<"Error in function algoDistEdgeApproxDominatingEdgesMessageBundling: \n";
        cout<<"Not enough memory to allocate the internal variables \n";
        exit(1);
    }
    PMessageBundle.resize(NumMessagesBundled*3, -1);//Initialize
    PCumulative.resize(numProcs+1, 0); //Only initialize the counter variable
    PSizeInfoMessages.resize(numProcs*3, 0);

    for (MilanInt i=0; i<numProcs; i++) // Changed by Fabio to be an integer, addresses needs to be integers!
        PCumulative[i+1]=PCumulative[i]+PCounter[i];
    //Reuse PCounter to keep track of how many messages were inserted:
    for (MilanInt i=0; i<numProcs; i++) // Changed by Fabio to be an integer, addresses needs to be integers!
        PCounter[i]=0;
    //Build the Message Bundle packet:
    for (MilanInt i=0; i<NumMessagesBundled; i++) { // Changed by Fabio to be an integer, addresses needs to be integers!
        myIndex = ( PCumulative[QOwner[i]] + PCounter[QOwner[i]] )*3;
        PMessageBundle[myIndex+0] = QLocalVtx[i];
        PMessageBundle[myIndex+1] = QGhostVtx[i];
        PMessageBundle[myIndex+2] = QMsgType[i];
        PCounter[QOwner[i]]++;
    }
    //Send the Bundled Messages: Use ISend
    vector<MPI_Request> SRequest; //Requests that are used for each send message
    vector<MPI_Status> SStatus;   //Status of sent messages, used in MPI_Wait
    MilanLongInt MessageIndex=0; //Pointer for current message
    try {
        SRequest.reserve(numProcs*2); //At most two messages per processor
        SStatus.reserve(numProcs*2);//At most two messages per processor
    } catch ( length_error ) {
        cout<<"Error in function algoDistEdgeApproxDominatingEdgesLinearSearchImmediateSend: \n";
        cout<<"Not enough memory to allocate the internal variables \n";
        exit(1);
    }
    MPI_Request myReq; //A sample request
    SRequest.resize(numProcs*2,myReq);
    MPI_Status myStat; //A sample status
    SStatus.resize(numProcs*2,myStat);
    //Send the Messages
    for (MilanInt i=0; i<numProcs; i++) { // Changed by Fabio to be an integer, addresses needs to be integers!
        if (i==myRank) //Do not send anything to yourself
            continue;
        //Send the Message with information about the size of next message:
        //Build the Message Packet:
        PSizeInfoMessages[i*3+0] = (PCumulative[i+1]-PCumulative[i])*3; // # of integers in the next message
        PSizeInfoMessages[i*3+1] = -1; //Dummy packet
        PSizeInfoMessages[i*3+2] = SIZEINFO;  //TYPE
        //Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
        cout<<"\n("<<myRank<<")Sending bundled message to process "<<i<<" size: "<<PSizeInfoMessages[i*3+0]<<endl;
        fflush(stdout);
#endif
        if ( PSizeInfoMessages[i*3+0] > 0 ) { //Send only if it is a nonempty packet
            MPI_Isend(&PSizeInfoMessages[i*3+0], 3, TypeMap<MilanLongInt>(), i, ComputeTag, comm, &SRequest[MessageIndex]);
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
            MPI_Isend(&PMessageBundle[PCumulative[i]*3], PSizeInfoMessages[i*3+0], TypeMap<MilanLongInt>(), i, BundleTag, comm, &SRequest[MessageIndex]);
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
    MilanInt OneMessageSize=0;
    MPI_Pack_size(3, TypeMap<MilanLongInt>(), comm, &OneMessageSize); //Size of one message packet
    //How many messages to send?
    //Potentially three kinds of messages will be sent/received:
    //Request, Success, Failure.
    //But only two will be sent from a given processor.
    //Substract the number of messages that have already been sent as bundled messages:
    MilanLongInt numMessagesToSend = numGhostEdges*2 - NumMessagesBundled;
    MilanInt     BufferSize = (OneMessageSize+MPI_BSEND_OVERHEAD)*numMessagesToSend;

    MilanLongInt *Buffer=0;
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
    ///////////////////////// END OF SEND BUNDLED MESSAGES //////////////////////////////////

    finishTime = MPI_Wtime();
    *ph1_time = finishTime-startTime; //Time taken for Phase-1
    *ph1_card = myCard ; //Cardinality at the end of Phase-1
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
                            msgInd++; msgActual++;
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

///Find the owner of a ghost node:
inline MilanInt findOwnerOfGhost(MilanLongInt vtxIndex, MilanLongInt *mVerDistance,
                                 MilanInt myRank, MilanInt numProcs) {
    //MilanLongInt Size = mVerDistance.size();
    MilanLongInt mStartInd = mVerDistance[myRank];
    MilanInt Start = 0;
    MilanInt End = numProcs;
    MilanInt Current = 0;

#if 0
    if ( vtxIndex < mStartInd )
    End = myRank;
  else
    Start = myRank;
#endif

    while ( Start <= End ) {
        Current = (End + Start)/2;
        //CASE-1:
        if ( mVerDistance[Current] == vtxIndex ) {
            while ( mVerDistance[Current+1] == vtxIndex ) {
                Current++;
                if ( Current == numProcs )
                    return (-1);
            }
            return (Current);
        }
        else { //CASE 2:
            if ( mVerDistance[Current] > vtxIndex )
                End = Current - 1;
            else //CASE 3:
                Start = Current + 1;
        }
    } //End of While()
    if ( Current == 0 )
        return (Current);
    else {
        if ( mVerDistance[Current] > vtxIndex )
            return (Current-1);
        else
            return (Current);
    } //End of else
    return (-1); //It should not reach here!
} //End of findOwnerOfGhost()

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
    for (int k = adj1; k < adj2; k++) {

        if ((edgeLocWeight[k] > heaviestEdgeWt) ||
            ((edgeLocWeight[k] == heaviestEdgeWt) && (w < verLocInd[k]))) {
            heaviestEdgeWt = edgeLocWeight[k];
            w = verLocInd[k];
        }
    } //End of for loop
    return w;
}

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
inline bool isAlreadyMatched(MilanLongInt k,
                                MilanLongInt* verLocInd,
                                MilanLongInt StartIndex,
                                MilanLongInt EndIndex,
                                vector <MilanLongInt> &GMate,
                                MilanLongInt* Mate,
                                map <MilanLongInt, MilanLongInt> &Ghost2LocalMap
                                ) {

    if ((verLocInd[k] < StartIndex) || (verLocInd[k] > EndIndex)) { //Is it a ghost vertex?
        if (GMate[Ghost2LocalMap[verLocInd[k]]] >= 0)// Already matched
            return true;
    } else { //A local vertex
        if (Mate[verLocInd[k] - StartIndex] >= 0) // Already matched
            return true;
    }

    return false;
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
                                              vector <MilanLongInt> &GMate,
                                              MilanLongInt* Mate,
                                              map <MilanLongInt, MilanLongInt> &Ghost2LocalMap)
{
    MilanInt w = -1;
    MilanReal heaviestEdgeWt = MilanRealMin; //Assign the smallest Value possible first LDBL_MIN
    for (k = adj1; k < adj2; k++) {
        if (isAlreadyMatched(k, verLocInd, StartIndex, EndIndex, GMate, Mate, Ghost2LocalMap)) continue;

        if ((edgeLocWeight[k] > heaviestEdgeWt) ||
            ((edgeLocWeight[k] == heaviestEdgeWt) && (w < verLocInd[k]))) {
            heaviestEdgeWt = edgeLocWeight[k];
            w = verLocInd[k];
        }
    } //End of for loop
    return w;
}
#endif

#endif