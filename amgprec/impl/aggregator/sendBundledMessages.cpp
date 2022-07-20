#include "MatchBoxPC.h"

void sendBundledMessages(MilanLongInt *numGhostEdges,
                                MilanInt *BufferSize,
                                MilanLongInt *Buffer,
                                vector<MilanLongInt> &PCumulative,
                                vector<MilanLongInt> &PMessageBundle,
                                vector<MilanLongInt> &PSizeInfoMessages,
                                MilanLongInt *PCounter,
                                MilanLongInt NumMessagesBundled,
                                MilanLongInt *msgActual,
                                MilanLongInt *msgInd,
                                MilanInt numProcs,
                                MilanInt myRank,
                                MPI_Comm comm,
                                vector<MilanLongInt> &QLocalVtx,
                                vector<MilanLongInt> &QGhostVtx,
                                vector<MilanLongInt> &QMsgType,
                                vector<MilanInt> &QOwner,
                                vector<MPI_Request> &SRequest,
                                vector<MPI_Status> &SStatus)
{

    MilanLongInt myIndex = 0, numMessagesToSend;
    MilanInt i = 0, OneMessageSize = 0;

#ifdef DEBUG_HANG_
    if (myRank == 0)
        cout << "\n(" << myRank << ") Send Bundles" << endl;
    fflush(stdout);
#endif

#pragma omp parallel private(i) default(shared) num_threads(NUM_THREAD)
    {
#pragma omp master
        {
// Data structures for Bundled Messages:
#pragma omp task depend(inout                                                       \
                        : PCumulative, PMessageBundle, PSizeInfoMessages) depend(in \
                                                                                 : NumMessagesBundled, numProcs)
            {try {
                PMessageBundle.reserve(NumMessagesBundled * 3); // Three integers per message
    PCumulative.reserve(numProcs + 1);                          // Similar to Row Pointer vector in CSR data structure
    PSizeInfoMessages.reserve(numProcs * 3);                    // Buffer to hold the Size info message packets
}
catch (length_error)
{
    cout << "Error in function algoDistEdgeApproxDominatingEdgesMessageBundling: \n";
    cout << "Not enough memory to allocate the internal variables \n";
    exit(1);
}
PMessageBundle.resize(NumMessagesBundled * 3, -1); // Initialize
PCumulative.resize(numProcs + 1, 0);               // Only initialize the counter variable
PSizeInfoMessages.resize(numProcs * 3, 0);
}

#pragma omp task depend(inout                    \
                        : PCumulative) depend(in \
                                              : PCounter)
{
    for (i = 0; i < numProcs; i++)
        PCumulative[i + 1] = PCumulative[i] + PCounter[i];
}

#pragma omp task depend(inout \
                        : PCounter)
{
    // Reuse PCounter to keep track of how many messages were inserted:
    for (MilanInt i = 0; i < numProcs; i++) // Changed by Fabio to be an integer, addresses needs to be integers!
        PCounter[i] = 0;
}

// Build the Message Bundle packet:
#pragma omp task depend(in                                                                                          \
                        : PCounter, QLocalVtx, QGhostVtx, QMsgType, QOwner, PMessageBundle, PCumulative) depend(out \
                                                                                                                : myIndex, PMessageBundle, PCounter)
{
    for (i = 0; i < NumMessagesBundled; i++)
    {
        myIndex = (PCumulative[QOwner[i]] + PCounter[QOwner[i]]) * 3;
        PMessageBundle[myIndex + 0] = QLocalVtx[i];
        PMessageBundle[myIndex + 1] = QGhostVtx[i];
        PMessageBundle[myIndex + 2] = QMsgType[i];
        PCounter[QOwner[i]]++;
    }
}

// Send the Bundled Messages: Use ISend
#pragma omp task depend(out \
                        : SRequest, SStatus)
{
    try
    {
        SRequest.reserve(numProcs * 2); // At most two messages per processor
        SStatus.reserve(numProcs * 2);  // At most two messages per processor
    }
    catch (length_error)
    {
        cout << "Error in function algoDistEdgeApproxDominatingEdgesLinearSearchImmediateSend: \n";
        cout << "Not enough memory to allocate the internal variables \n";
        exit(1);
    }
}

// Send the Messages
#pragma omp task depend(inout                                                  \
                        : SRequest, PSizeInfoMessages, PCumulative) depend(out \
                                                                           : *msgActual, *msgInd)
{
    for (i = 0; i < numProcs; i++)
    {                    // Changed by Fabio to be an integer, addresses needs to be integers!
        if (i == myRank) // Do not send anything to yourself
            continue;
        // Send the Message with information about the size of next message:
        // Build the Message Packet:
        PSizeInfoMessages[i * 3 + 0] = (PCumulative[i + 1] - PCumulative[i]) * 3; // # of integers in the next message
        PSizeInfoMessages[i * 3 + 1] = -1;                                        // Dummy packet
        PSizeInfoMessages[i * 3 + 2] = SIZEINFO;                                  // TYPE
                                                                                  // Send a Request (Asynchronous)
#ifdef PRINT_DEBUG_INFO_
        cout << "\n(" << myRank << ")Sending bundled message to process " << i << " size: " << PSizeInfoMessages[i * 3 + 0] << endl;
        fflush(stdout);
#endif
        if (PSizeInfoMessages[i * 3 + 0] > 0)
        { // Send only if it is a nonempty packet
            MPI_Isend(&PSizeInfoMessages[i * 3 + 0], 3, TypeMap<MilanLongInt>(), i, ComputeTag, comm,
                      &SRequest[(*msgInd)]);
            (*msgActual)++;
            (*msgInd)++;
            // Now Send the message with the data packet:
#ifdef PRINT_DEBUG_INFO_
            cout << "\n(" << myRank << ")SendiFFng Bundle to : " << i << endl;
            for (k = (PCumulative[i] * 3); k < (PCumulative[i] * 3 + PSizeInfoMessages[i * 3 + 0]); k++)
                cout << PMessageBundle[k] << ",";
            cout << endl;
            fflush(stdout);
#endif
            MPI_Isend(&PMessageBundle[PCumulative[i] * 3], PSizeInfoMessages[i * 3 + 0],
                      TypeMap<MilanLongInt>(), i, BundleTag, comm, &SRequest[(*msgInd)]);
            (*msgInd)++;
        } // End of if size > 0
    }
}

#pragma omp task depend(inout \
                        : PCumulative, QLocalVtx, QGhostVtx, QMsgType, QOwner)
{

    // Free up temporary memory:
    PCumulative.clear();
    QLocalVtx.clear();
    QGhostVtx.clear();
    QMsgType.clear();
    QOwner.clear();
}

#pragma omp task depend(inout : OneMessageSize, *BufferSize) depend(out : numMessagesToSend) depend(in : *numGhostEdges)
{

#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << ")Number of Ghost edges = " << *numGhostEdges;
    cout << "\n(" << myRank << ")Total number of potential message X 2 = " << *numGhostEdges * 2;
    cout << "\n(" << myRank << ")Number messages already sent in bundles = " << NumMessagesBundled;
    if (*numGhostEdges > 0)
    {
        cout << "\n(" << myRank << ")Percentage of total = " << ((double)NumMessagesBundled / (double)(*numGhostEdges * 2)) * 100.0 << "% \n";
    }
    fflush(stdout);
#endif

    // Allocate memory for MPI Send messages:
    /* WILL COME BACK HERE - NO NEED TO STORE ALL THIS MEMORY !! */
    OneMessageSize = 0;
    MPI_Pack_size(3, TypeMap<MilanLongInt>(), comm, &OneMessageSize); // Size of one message packet
    // How many messages to send?
    // Potentially three kinds of messages will be sent/received:
    // Request, Success, Failure.
    // But only two will be sent from a given processor.
    // Substract the number of messages that have already been sent as bundled messages:
    numMessagesToSend = (*numGhostEdges) * 2 - NumMessagesBundled;
    *BufferSize = (OneMessageSize + MPI_BSEND_OVERHEAD) * numMessagesToSend;
}

#pragma omp task depend(out : Buffer) depend(in : *BufferSize)
{
    Buffer = 0;
#ifdef PRINT_DEBUG_INFO_
    cout << "\n(" << myRank << ")Size of One Message from PACK= " << OneMessageSize;
    cout << "\n(" << myRank << ")Size of Message overhead = " << MPI_BSEND_OVERHEAD;
    cout << "\n(" << myRank << ")Number of Ghost edges = " << *numGhostEdges;
    cout << "\n(" << myRank << ")Number of remaining message = " << numMessagesToSend;
    cout << "\n(" << myRank << ")BufferSize = " << (*BufferSize);
    cout << "\n(" << myRank << ")Attaching Buffer on.. ";
    fflush(stdout);
#endif
    if ((*BufferSize) > 0)
    {
        Buffer = (MilanLongInt *)malloc((*BufferSize)); // Allocate memory
        if (Buffer == 0)
        {
            cout << "Error in function algoDistEdgeApproxDominatingEdgesLinearSearch: \n";
            cout << "Not enough memory to allocate for send buffer on process " << myRank << "\n";
            exit(1);
        }
        MPI_Buffer_attach(Buffer, *BufferSize); // Attach the Buffer
    }
}
}
}
}