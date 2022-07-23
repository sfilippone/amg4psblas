#include "MatchBoxPC.h"

void queuesTransfer(staticQueue &U,
                    staticQueue &privateU,
                    vector<MilanLongInt> &QLocalVtx,
                    vector<MilanLongInt> &QGhostVtx,
                    vector<MilanLongInt> &QMsgType,
                    vector<MilanInt> &QOwner,
                    vector<MilanLongInt> &privateQLocalVtx,
                    vector<MilanLongInt> &privateQGhostVtx,
                    vector<MilanLongInt> &privateQMsgType,
                    vector<MilanInt> &privateQOwner)
{

#pragma omp critical(U)
    {
        while (!privateU.empty())
            U.push_back(privateU.pop_back());
    }

#pragma omp critical(sendMessageTransfer)
    {

        QLocalVtx.insert(QLocalVtx.end(), privateQLocalVtx.begin(), privateQLocalVtx.end());
        QGhostVtx.insert(QGhostVtx.end(), privateQGhostVtx.begin(), privateQGhostVtx.end());
        QMsgType.insert(QMsgType.end(), privateQMsgType.begin(), privateQMsgType.end());
        QOwner.insert(QOwner.end(), privateQOwner.begin(), privateQOwner.end());
    }

    privateQLocalVtx.clear();
    privateQGhostVtx.clear();
    privateQMsgType.clear();
    privateQOwner.clear();
}