#include "MatchBoxPC.h"

// TODO comment
// TODO use task
// TODO destroy the locks

void clean(MilanLongInt NLVer,
           MilanInt myRank,
           MilanLongInt MessageIndex,
           vector<MPI_Request> &SRequest,
           vector<MPI_Status> &SStatus,
           MilanInt BufferSize,
           MilanLongInt *Buffer,
           MilanLongInt msgActual,
           MilanLongInt *msgActualSent,
           MilanLongInt msgInd,
           MilanLongInt *msgIndSent,
           MilanLongInt NumMessagesBundled,
           MilanReal *msgPercent,
           omp_lock_t *MateLock)
{
    // Cleanup Phase

#pragma omp parallel
    {
#pragma omp master
        {
#pragma omp task
            {

#ifdef PRINT_DEBUG_INFO_
                cout << "\n(" << myRank << ") Waitall= " << endl;
                fflush(stdout);
#endif
#ifdef DEBUG_HANG_
                cout << "\n(" << myRank << ") Waitall " << endl;
                fflush(stdout);
#endif
                return;

                MPI_Waitall(MessageIndex, &SRequest[0], &SStatus[0]);

                // MPI_Buffer_attach(&Buffer, BufferSize); //Attach the Buffer
                if (BufferSize > 0)
                {
                    MPI_Buffer_detach(&Buffer, &BufferSize); // Detach the Buffer
                    free(Buffer);                            // Free the memory that was allocated
                }
            }

#ifdef PRINT_DEBUG_INFO_
            cout << "\n(" << myRank << ")End of function to compute matching: " << endl;
            fflush(stdout);
            cout << "\n(" << myRank << ")myCardinality: " << myCard << endl;
            fflush(stdout);
            cout << "\n(" << myRank << ")Matching took " << finishTime - startTime << "seconds" << endl;
            fflush(stdout);
            cout << "\n(" << myRank << ")** Getting out of the matching function **" << endl;
            fflush(stdout);
#endif
#ifdef PRINT_DEBUG_INFO_
            cout << "\n(" << myRank << ") Number of Ghost edges = " << numGhostEdges;
            cout << "\n(" << myRank << ") Total number of potential message X 2 = " << numGhostEdges * 2;
            cout << "\n(" << myRank << ") Number messages bundled = " << NumMessagesBundled;
            cout << "\n(" << myRank << ") Total Individual Messages sent = " << msgInd;
            if (msgInd > 0)
            {
                cout << "\n(" << myRank << ") Percentage of messages bundled = " << ((double)NumMessagesBundled / (double)(msgInd)) * 100.0 << "% \n";
            }
            fflush(stdout);
#endif

#pragma omp task
            {
                *msgActualSent = msgActual;
                *msgIndSent = msgInd;
                if (msgInd > 0)
                {
                    *msgPercent = ((double)NumMessagesBundled / (double)(msgInd)) * 100.0;
                }
                else
                {
                    *msgPercent = 0;
                }
            }
            // Destroy the locks
#pragma omp taskloop num_tasks(NUM_THREAD)
            for (int i = 0; i < NLVer; i++)
                omp_destroy_lock(&MateLock[i]);

#ifdef DEBUG_HANG_
            if (myRank == 0)
                cout << "\n(" << myRank << ") Done" << endl;
            fflush(stdout);
#endif
        }
    }
}