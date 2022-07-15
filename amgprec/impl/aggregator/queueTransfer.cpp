#include "MatchBoxPC.h"

void queuesTransfer(staticQueue &U,
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
            U.push_back(privateU.pop_back());
    }

#pragma omp critical(privateMsg)
    {
        while (!privateQLocalVtx.empty())
        {
            QLocalVtx.push_back(privateQLocalVtx.pop_back());
            QGhostVtx.push_back(privateQGhostVtx.pop_back());
            QMsgType.push_back(privateQMsgType.pop_back());
            QOwner.push_back(privateQOwner.pop_back());
        }
    }
}