#include "MatchBoxPC.h"
#include <stdio.h>
#include <iostream>
#include <map>
#include <vector>
#include "primitiveDataTypeDefinitions.h"
#include "dataStrStaticQueue.h"
#include "omp.h"

inline void processMessages(int error_codeC,
                            MilanInt numProcs,
                            MilanInt myRank,
                            int ComputeTag,
                            int BundleTag,
                            MPI_Comm comm,
                            vector<MilanLongInt> &Message,
                            char *error_message,
                            int message_length,
                            vector<MilanLongInt> &ReceiveBuffer,
                            MilanLongInt *BundleSizePtr)
{

    MilanInt Sender;
    MPI_Status computeStatus;
    MilanLongInt bundleSize = *BundleSizePtr;

#ifdef PRINT_DEBUG_INFO_
    cout
        << "\n(" << myRank << "=========================************===============================" << endl;
    fflush(stdout);
    fflush(stdout);
#endif
#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << ")About to begin Message processing phase ... S=" << S << endl;
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

    *BundleSizePtr = bundleSize;
    return;
}