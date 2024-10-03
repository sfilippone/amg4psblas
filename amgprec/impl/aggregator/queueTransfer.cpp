#include "MatchBoxPC.h"

void queuesTransfer(vector<MilanLongInt> &U,
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

#pragma omp critical(U)
    {
        U.insert(U.end(), privateU.begin(), privateU.end());
    }

    privateU.clear();

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
