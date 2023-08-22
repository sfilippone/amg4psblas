#include "MatchBoxPC.h"

void PARALLEL_PROCESS_EXPOSED_VERTEX_B(MilanLongInt NLVer,
                                       MilanLongInt *candidateMate,
                                       MilanLongInt *verLocInd,
                                       MilanLongInt *verLocPtr,
                                       MilanLongInt StartIndex,
                                       MilanLongInt EndIndex,
                                       MilanLongInt *Mate,
                                       vector<MilanLongInt> &GMate,
                                       map<MilanLongInt, MilanLongInt> &Ghost2LocalMap,
                                       MilanReal *edgeLocWeight,
                                       MilanLongInt *myCard,
                                       MilanLongInt *msgInd,
                                       MilanLongInt *NumMessagesBundled,
                                       MilanLongInt *S,
                                       MilanLongInt *verDistance,
                                       MilanLongInt *PCounter,
                                       vector<MilanLongInt> &Counter,
                                       MilanInt myRank,
                                       MilanInt numProcs,
                                       vector<MilanLongInt> &U,
                                       vector<MilanLongInt> &privateU,
                                       vector<MilanLongInt> &QLocalVtx,
                                       vector<MilanLongInt> &QGhostVtx,
                                       vector<MilanLongInt> &QMsgType,
                                       vector<MilanInt> &QOwner,
                                       vector<MilanLongInt> &privateQLocalVtx,
                                       vector<MilanLongInt> &privateQGhostVtx,
                                       vector<MilanLongInt> &privateQMsgType,
                                       vector<MilanInt> &privateQOwner)
{ 

    MilanLongInt v = -1, k = -1, w = -1, adj11 = 0, adj12 = 0, k1 = 0;
    MilanInt ghostOwner = 0, option, igw;

#pragma omp parallel private(option, k, w, v, k1, adj11, adj12, ghostOwner)    \
    firstprivate(privateU, StartIndex, EndIndex, privateQLocalVtx, privateQGhostVtx, privateQMsgType, privateQOwner) \
     default(shared)   num_threads(NUM_THREAD)

    {
#pragma omp for reduction(+                             \
                          : PCounter[:numProcs], myCard \
                          [:1], msgInd                  \
                          [:1], NumMessagesBundled      \
                          [:1]) \
        schedule(static)
        for (v = 0; v < NLVer; v++) {
            option = -1;
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

#pragma omp critical(processExposed)
                {
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
		      (*myCard)++;
		      if ((w < StartIndex) || (w > EndIndex)) { // w is a ghost vertex
			option = 2;
			if (candidateMate[NLVer + Ghost2LocalMap[w]] == v + StartIndex) {
			  option = 1;
			  Mate[v] = w;
			  GMate[Ghost2LocalMap[w]] = v + StartIndex; // w is a Ghost
			  
			} // End of if CandidateMate[w] = v

		      } // End of if a Ghost Vertex
		      else  { // w is a local vertex
			
			if (candidateMate[w - StartIndex] == (v + StartIndex)) {
			  option = 3;
			  Mate[v] = w;                           // v is local
			  Mate[w - StartIndex] = v + StartIndex; // w is local
			  
#ifdef PRINT_DEBUG_INFO_
			  cout << "\n(" << myRank << ")MATCH: (" << v + StartIndex << "," << w << ") ";
			  fflush(stdout);
#endif
			  
			} // End of if ( candidateMate[w-StartIndex] == (v+StartIndex) )
		      }     // End of Else
		      
                    } // End of second if
		    
                } // End critical processExposed
		
            } // End of if(w >=0)
            else  {
	      // This piece of code is executed a really small amount of times
	      adj11 = verLocPtr[v];
	      adj12 = verLocPtr[v + 1];
	      for (k1 = adj11; k1 < adj12; k1++) {
		w = verLocInd[k1];
		if ((w < StartIndex) || (w > EndIndex)) { // A ghost

#ifdef PRINT_DEBUG_INFO_
		  cout << "\n(" << myRank << ")Sending a failure message: ";
		  cout << "\n(" << myRank << ")Ghost is " << w << " Owner is: " << findOwnerOfGhost(w, verDistance, myRank, numProcs);
		  fflush(stdout);
#endif
		  (*msgInd)++;
		  (*NumMessagesBundled)++;
		  ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
		  // assert(ghostOwner != -1);
		  // assert(ghostOwner != myRank);
		  PCounter[ghostOwner]++;
		  
		  privateQLocalVtx.push_back(v + StartIndex);
		  privateQGhostVtx.push_back(w);
		  privateQMsgType.push_back(FAILURE);
		  privateQOwner.push_back(ghostOwner);
		  
		} // End of if(GHOST)
	      }     // End of for loop
            }
            // End:   PARALLEL_PROCESS_EXPOSED_VERTEX_B(v)
	    
            switch (option)
            {
            case -1:
                break;
            case 1:
                privateU.push_back(v + StartIndex);
                privateU.push_back(w);

#ifdef PRINT_DEBUG_INFO_
                cout << "\n(" << myRank << ")MATCH: (" << v + StartIndex << "," << w << ")";
                fflush(stdout);
#endif

                //  Decrement the counter:
                PROCESS_CROSS_EDGE(&Counter[Ghost2LocalMap[w]], S);
            case 2:
#ifdef PRINT_DEBUG_INFO_
                cout << "\n(" << myRank << ")Sending a request message (291):";
                cout << "\n(" << myRank << ")Local is: " << v + StartIndex << " Ghost is " << w << " Owner is: " << findOwnerOfGhost(w, verDistance, myRank, numProcs) << endl;
                fflush(stdout);
#endif
                (*msgInd)++;
                (*NumMessagesBundled)++;
                ghostOwner = findOwnerOfGhost(w, verDistance, myRank, numProcs);
                // assert(ghostOwner != -1);
                // assert(ghostOwner != myRank);
                PCounter[ghostOwner]++;

                privateQLocalVtx.push_back(v + StartIndex);
                privateQGhostVtx.push_back(w);
                privateQMsgType.push_back(REQUEST);
                privateQOwner.push_back(ghostOwner);
                break;
            case 3:
            default:
                privateU.push_back(v + StartIndex);
                privateU.push_back(w);
                break;
            }

        } // End of for ( v=0; v < NLVer; v++ )

        queuesTransfer(U, privateU, QLocalVtx,
                       QGhostVtx,
                       QMsgType, QOwner, privateQLocalVtx,
                       privateQGhostVtx,
                       privateQMsgType,
                       privateQOwner);

    } // End of parallel region
}
