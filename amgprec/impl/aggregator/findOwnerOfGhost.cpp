#include "MatchBoxPC.h"

/// Find the owner of a ghost node:
MilanInt findOwnerOfGhost(MilanLongInt vtxIndex, MilanLongInt *mVerDistance,
                          MilanInt myRank, MilanInt numProcs)
{

  MilanLongInt mStartInd = mVerDistance[myRank];
  MilanInt Start = 0;
  MilanInt End = numProcs;
  MilanInt Current = 0;

  while (Start <= End)
  {
    Current = (End + Start) / 2;
    // CASE-1:
    if (mVerDistance[Current] == vtxIndex) return Current;
    else // CASE 2:
      if (mVerDistance[Current] > vtxIndex)
        End = Current - 1;
      else // CASE 3:
        Start = Current + 1;
  } // End of While()

 if (mVerDistance[Current] > vtxIndex)
      return (Current - 1);

  return Current;
} // End of findOwnerOfGhost()
