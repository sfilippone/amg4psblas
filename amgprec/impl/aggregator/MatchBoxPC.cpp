// ***********************************************************************
//
//        MatchboxP: A C++ library for approximate weighted matching
//               Mahantesh Halappanavar (hala@pnnl.gov)
//               Pacific Northwest National Laboratory
//
// ***********************************************************************
//
//       Copyright (2021) Battelle Memorial Institute
//                      All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "MatchBoxPC.h"
#ifdef __cplusplus
extern "C" {
#endif


void dMatchBoxPC(MilanLongInt NLVer, MilanLongInt NLEdge,
		MilanLongInt* verLocPtr, MilanLongInt* verLocInd, MilanReal* edgeLocWeight,
		MilanLongInt* verDistance,
		MilanLongInt* Mate,
		MilanInt myRank, MilanInt numProcs, MilanInt icomm,
		MilanLongInt* msgIndSent, MilanLongInt* msgActualSent, MilanReal* msgPercent,
		MilanReal* ph0_time, MilanReal* ph1_time, MilanReal* ph2_time,
		MilanLongInt* ph1_card, MilanLongInt* ph2_card ) {
  MPI_Comm C_comm=MPI_Comm_f2c(icomm);
#ifdef DEBUG
  fprintf(stderr,"MatchBoxPC: rank %d nlver %ld nledge %ld [ %ld %ld ]\n",
	  myRank,NLVer, NLEdge,verDistance[0],verDistance[1]);
#endif
  dalgoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMateC(NLVer, NLEdge,
							   verLocPtr, verLocInd, edgeLocWeight,
							   verDistance,  Mate,
							   myRank, numProcs, C_comm,
							   msgIndSent, msgActualSent, msgPercent,
							   ph0_time, ph1_time, ph2_time,
							   ph1_card, ph2_card );
}

void sMatchBoxPC(MilanLongInt NLVer, MilanLongInt NLEdge,
		MilanLongInt* verLocPtr, MilanLongInt* verLocInd, MilanFloat* edgeLocWeight,
		MilanLongInt* verDistance,
		MilanLongInt* Mate,
		MilanInt myRank, MilanInt numProcs, MilanInt icomm,
		MilanLongInt* msgIndSent, MilanLongInt* msgActualSent, MilanReal* msgPercent,
		MilanReal* ph0_time, MilanReal* ph1_time, MilanReal* ph2_time,
		MilanLongInt* ph1_card, MilanLongInt* ph2_card ) {
  MPI_Comm C_comm=MPI_Comm_f2c(icomm);
#ifdef DEBUG
  fprintf(stderr,"MatchBoxPC: rank %d nlver %ld nledge %ld [ %ld %ld ]\n",
	  myRank,NLVer, NLEdge,verDistance[0],verDistance[1]);
#endif
  salgoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMateC(NLVer, NLEdge,
							   verLocPtr, verLocInd, edgeLocWeight,
							   verDistance,  Mate,
							   myRank, numProcs, C_comm,
							   msgIndSent, msgActualSent, msgPercent,
							   ph0_time, ph1_time, ph2_time,
							   ph1_card, ph2_card );
}

#ifdef __cplusplus
}
#endif
