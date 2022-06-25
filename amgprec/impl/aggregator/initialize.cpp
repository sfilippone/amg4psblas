#include "MatchBoxPC.h"
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <map>
#include <vector>
#include "primitiveDataTypeDefinitions.h"
#include "dataStrStaticQueue.h"
#include "omp.h"

inline void initialize(MilanLongInt NLVer, MilanLongInt NLEdge,
                        MilanLongInt StartIndex, MilanLongInt EndIndex,
                        MilanLongInt* numGhostEdgesPtr,
                        MilanLongInt* numGhostVerticesPtr,
                        MilanLongInt* insertMePtr,
                        MilanLongInt* verLocInd,
                        MilanLongInt* verLocPtr,
                        omp_lock_t* MateLock,
                        map <MilanLongInt, MilanLongInt> &Ghost2LocalMap,
                        vector <MilanLongInt>& Counter,
                        vector <MilanLongInt>& verGhostPtr,
                        vector <MilanLongInt>& verGhostInd,
                        vector <MilanLongInt>& tempCounter,
                        vector <MilanLongInt>& GMate,
                        vector<MilanLongInt>& Message,
                        vector<MilanLongInt>& QLocalVtx,
                        vector<MilanLongInt>& QGhostVtx,
                        vector<MilanLongInt>& QMsgType,
                        vector<MilanInt>& QOwner,
                        MilanLongInt* candidateMate,
                        staticQueue& U
                        )
{

    MilanLongInt insertMe = 0, numGhostEdges = 0, numGhostVertices = 0;
    MilanLongInt adj1, adj2;
    int i, v, k, w;

    
    // index that starts with zero to |Vg|  - 1
    map<MilanLongInt, MilanLongInt>::iterator storedAlready;

#pragma omp parallel private(insertMe, k, w, v, adj1, adj2) firstprivate(StartIndex, EndIndex) default(shared) num_threads(4)
    {
        
        //Initialize the locks
        //TODO this can be executed as task in parallel with other unparallelizable tasks
        //TODO destroy the locks
#pragma omp for schedule(static)
        for(i = 0; i < NLVer; i++)
            omp_init_lock(&MateLock[i]);


#ifdef TIME_TRACKER
    double Ghost2LocalInitialization = MPI_Wtime();
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


            /*
             * Not parallelizable
             */

            for (i = 0; i < numGhostVertices; i++) { //O(|Ghost Vertices|)
                verGhostPtr[i + 1] = verGhostPtr[i] + Counter[i];
#ifdef PRINT_DEBUG_INFO_
                cout<<verGhostPtr[i]<<"\t"; fflush(stdout);
#endif
            }
        } // End of single region

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

#pragma omp for nowait schedule(static)
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
    
    }

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
            //message_type = 0;
            //NumMessagesBundled = 0;
            //ghostOwner = 0;
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

        } // end of single region

    *numGhostEdgesPtr = numGhostEdges;
    *numGhostVerticesPtr = numGhostVertices;  
    *insertMePtr = insertMe; 
}
