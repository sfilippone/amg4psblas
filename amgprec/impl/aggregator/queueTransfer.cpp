#include "MatchBoxPC.h"
#include <stdio.h>
#include <iostream>
#include <map>
#include <vector>
#include "primitiveDataTypeDefinitions.h"
#include "dataStrStaticQueue.h"
#include "omp.h"

inline void queuesTransfer(staticQueue &U,
                           staticQueue &privateU,
                           vector<MilanLongInt> &QLocalVtx,
                           vector<MilanLongInt> &QGhostVtx,
                           vector<MilanLongInt> &QMsgType,
                           vector<MilanInt> &QOwner,
                           staticQueue &privateQLocalVtx,
                           staticQueue &privateQGhostVtx,
                           staticQueue &privateQMsgType,
                           staticQueue &privateQOwner)
{

#pragma omp critical(U)
    {
        while (!privateU.empty())
            U.push_back(privateU.pop_front());
    }

#pragma omp critical(privateMsg)
    {
        while (!privateQLocalVtx.empty())
        {
            QLocalVtx.push_back(privateQLocalVtx.pop_front());
            QGhostVtx.push_back(privateQGhostVtx.pop_front());
            QMsgType.push_back(privateQMsgType.pop_front());
            QOwner.push_back(privateQOwner.pop_front());
        }
    }
}