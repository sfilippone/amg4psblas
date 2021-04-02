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

#ifndef _primitiveDataType_Definition_
#define _primitiveDataType_Definition_

#include "preProcessorDirectives.h"

using namespace std;

//Comment out these if you do not need 64 bits.
//#ifndef BIT64
//	#define BIT64
//#endif

//Regular integer:
#ifndef INTEGER_H
#define INTEGER_H
	typedef int MilanInt;
//	typedef MPI_INT MilanMpiInt;
#endif

//Regular long Integer:
#ifndef LONG_INT_H
#define LONG_INT_H
	#ifdef BIT64
	typedef int64_t MilanLongInt;
//	typedef MPI_LONG MilanMpiLongInt;
	#else
	typedef int MilanLongInt;
//	typedef MPI_INT MilanMpiLongInt;
	#endif
#endif

//Regular boolean
#ifndef BOOL_H
#define BOOL_H
	typedef bool MilanBool;
#endif

//Regular double and the Absolute Function:
#ifndef REAL_H
#define REAL_H
	typedef double MilanReal;
	//typedef MPI_DOUBLE MilanMpiReal;
	inline MilanReal MilanAbs(MilanReal value)
	{
	  return fabs(value);
	}
#endif

//Regular double and the Absolute Function:
#ifndef FLOAT_H
#define FLOAT_H
	typedef float MilanFloat;
	//typedef MPI_FLOAT MilanMpiFloat;
	inline MilanFloat MilanAbsFloat(MilanFloat value)
	{
	  return fabs(value);
	}
#endif

//// Define the limits:
#ifndef LIMITS_H
#define LIMITS_H

//Integer Maximum and Minimum:
#define MilanIntMax INT_MAX
#define MilanIntMin INT_MIN

#ifdef BIT64
	#define MilanLongIntMax LONG_MAX
	#define MilanLongIntMin -LONG_MAX
#else
	#define MilanLongIntMax INT_MAX
	#define MilanLongIntMin -INT_MAX
#endif

//Double Maximum and Minimum:
//Note: You can alternative use INFINITY defined in math.h
//It has been my experience that this is not very portable.
//Therefore I have adopted for LDBL_MAX and LDBL_MIN as +/- infinity.

//Largest positive number: LDBL_MAX = +infinity
//Smallest positive number: LDBL_MIN
//Smallest negative number: -LDBL_MAX = -infinity
//Largest negative number: -LDBL_MIN  (just next to zero on the other side?)

// +INFINITY
const double PLUS_INFINITY = numeric_limits<double>::infinity();
const float FPLUS_INFINITY = numeric_limits<float>::infinity();
//if(numeric_limits<float>::has_infinity)
// PLUS_INFINITY=numeric_limits<float>::infinity();
//else cerr<<"infinity for float isn�t supported";

const double MINUS_INFINITY = -PLUS_INFINITY;
const float FMINUS_INFINITY = -FPLUS_INFINITY;


//#define MilanRealMax LDBL_MAX
#define MilanRealMax PLUS_INFINITY
#define MilanFloatMax FPLUS_INFINITY

// -INFINITY
//Instead of assigning smallest possible positive number, assign smallest negative number
//although we only consider postive weights, just for correctness of understand.
//#define MilanRealMin -LDBL_MAX
//#define MilanRealMin LDBL_MIN
#define MilanRealMin MINUS_INFINITY
#define MilanFloatMin FMINUS_INFINITY

//const double PLUS_INFINITY = LDBL_MAX;   //deprecated



#endif

#endif
