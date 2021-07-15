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

/*
 Feature: Message Aggregation:
 Request messages from the Initialization phase are aggregated and only
 one message per processor is sent.
 Data structures for message aggregation are similar to the data structures
 for storing compressed matrices/graphs - Pointers + actual messages.
 Assumption: processor indices are numbered from zero to P-1.
 */

/* Special feature: Mate is proportional to the size of local number of verices */

#ifndef _matchboxpC_H_
#define _matchboxpC_H_
//Turn on a lot of debugging information with this switch:
//#define PRINT_DEBUG_INFO_
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <map>
#include <vector>
// #include "matchboxp.h"
#include "primitiveDataTypeDefinitions.h"
#include "dataStrStaticQueue.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(SERIAL_MPI)
  
#define MilanMpiLongInt  MPI_LONG_LONG

#ifndef _primitiveDataType_Definition_
#define _primitiveDataType_Definition_
    //Regular integer:
    #ifndef INTEGER_H
    #define INTEGER_H
        typedef int32_t MilanInt;
    #endif

    //Regular long integer:
    #ifndef LONG_INT_H
    #define LONG_INT_H
        #ifdef BIT64
            typedef int64_t MilanLongInt;
            typedef MPI_LONG MilanMpiLongInt;
        #else
            typedef int32_t MilanLongInt;
            typedef MPI_INT MilanMpiLongInt;
        #endif
    #endif

    //Regular boolean
    #ifndef BOOL_H
    #define BOOL_H
        typedef bool MilanBool;
    #endif

    //Regular double and absolute value computation:
    #ifndef REAL_H
    #define REAL_H
        typedef double MilanReal;
        typedef MPI_DOUBLE MilanMpiReal;
        inline MilanReal MilanAbs(MilanReal value)
        {
            return fabs(value);
        }
    #endif

    //Regular float and absolute value computation:
    #ifndef FLOAT_H
    #define FLOAT_H
        typedef float MilanFloat;
        typedef MPI_FLOAT MilanMpiFloat;
        inline MilanFloat MilanAbsFloat(MilanFloat value)
        {
            return fabs(value);
        }
    #endif

    //// Define the limits:
    #ifndef LIMITS_H
    #define LIMITS_H
    //Integer Maximum and Minimum:
  //      #define MilanIntMax INT_MAX
  //    #define MilanIntMin INT_MIN
        #define MilanIntMax INT32_MAX
        #define MilanIntMin INT32_MIN

        #ifdef BIT64
            #define MilanLongIntMax INT64_MAX
            #define MilanLongIntMin -INT64_MAX
        #else
            #define MilanLongIntMax INT32_MAX
            #define MilanLongIntMin -INT32_MAX
        #endif

    #endif

    // +INFINITY
    const double PLUS_INFINITY = numeric_limits<int>::infinity();
    const double MINUS_INFINITY = -PLUS_INFINITY;
    //#define MilanRealMax LDBL_MAX
    #define MilanRealMax PLUS_INFINITY
    #define MilanRealMin MINUS_INFINITY
#endif

//Function of find the owner of a ghost vertex using binary search:
inline MilanInt findOwnerOfGhost(MilanLongInt vtxIndex, MilanLongInt *mVerDistance,
                                     MilanInt myRank, MilanInt numProcs);

  void dalgoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMateC
(
 MilanLongInt NLVer, MilanLongInt NLEdge,
 MilanLongInt* verLocPtr, MilanLongInt* verLocInd, MilanReal* edgeLocWeight,
 MilanLongInt* verDistance,
 MilanLongInt* Mate,
 MilanInt myRank, MilanInt numProcs, MPI_Comm comm,
 MilanLongInt* msgIndSent, MilanLongInt* msgActualSent, MilanReal* msgPercent,
 MilanReal* ph0_time, MilanReal* ph1_time, MilanReal* ph2_time,
 MilanLongInt* ph1_card, MilanLongInt* ph2_card );

 void salgoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMateC
(
MilanLongInt NLVer, MilanLongInt NLEdge,
MilanLongInt* verLocPtr, MilanLongInt* verLocInd, MilanFloat* edgeLocWeight,
MilanLongInt* verDistance,
MilanLongInt* Mate,
MilanInt myRank, MilanInt numProcs, MPI_Comm comm,
MilanLongInt* msgIndSent, MilanLongInt* msgActualSent, MilanReal* msgPercent,
MilanReal* ph0_time, MilanReal* ph1_time, MilanReal* ph2_time,
MilanLongInt* ph1_card, MilanLongInt* ph2_card );

void dMatchBoxPC(MilanLongInt NLVer, MilanLongInt NLEdge,
		MilanLongInt* verLocPtr, MilanLongInt* verLocInd, MilanReal* edgeLocWeight,
		MilanLongInt* verDistance,
		MilanLongInt* Mate,
		MilanInt myRank, MilanInt numProcs, MilanInt icomm,
		MilanLongInt* msgIndSent, MilanLongInt* msgActualSent, MilanReal* msgPercent,
		MilanReal* ph0_time, MilanReal* ph1_time, MilanReal* ph2_time,
		MilanLongInt* ph1_card, MilanLongInt* ph2_card );

void sMatchBoxPC(MilanLongInt NLVer, MilanLongInt NLEdge,
		MilanLongInt* verLocPtr, MilanLongInt* verLocInd, MilanFloat* edgeLocWeight,
		MilanLongInt* verDistance,
		MilanLongInt* Mate,
		MilanInt myRank, MilanInt numProcs, MilanInt icomm,
		MilanLongInt* msgIndSent, MilanLongInt* msgActualSent, MilanReal* msgPercent,
		MilanReal* ph0_time, MilanReal* ph1_time, MilanReal* ph2_time,
		MilanLongInt* ph1_card, MilanLongInt* ph2_card );

#endif
#ifdef __cplusplus
}
#endif
#endif