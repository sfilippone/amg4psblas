#include "MatchBoxPC.h"
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include "primitiveDataTypeDefinitions.h"
#include "dataStrStaticQueue.h"

///Find the owner of a ghost node:
inline MilanInt findOwnerOfGhost(MilanLongInt vtxIndex, MilanLongInt *mVerDistance,
                                     MilanInt myRank, MilanInt numProcs) {
  //MilanLongInt Size = mVerDistance.size();
  MilanLongInt mStartInd = mVerDistance[myRank];
  MilanInt Start = 0;
  MilanInt End = numProcs;
  MilanInt Current = 0;

#if 0
  if ( vtxIndex < mStartInd )
    End = myRank;
  else
    Start = myRank;
#endif

  while ( Start <= End ) {
    Current = (End + Start)/2;
    //CASE-1:
    if ( mVerDistance[Current] == vtxIndex ) {
      while ( mVerDistance[Current+1] == vtxIndex ) {
	Current++;
	if ( Current == numProcs )
	  return (-1);
      }
      return (Current);
    }
    else { //CASE 2:
      if ( mVerDistance[Current] > vtxIndex )
	End = Current - 1;
      else //CASE 3:
	Start = Current + 1;
    }
  } //End of While()
  if ( Current == 0 )
    return (Current);
  else {
    if ( mVerDistance[Current] > vtxIndex )
      return (Current-1);
    else
      return (Current);
  } //End of else
  return (-1); //It should not reach here!
} //End of findOwnerOfGhost()
