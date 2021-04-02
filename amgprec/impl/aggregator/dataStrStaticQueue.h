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

#ifndef _static_Queue_
#define _static_Queue_

#include "primitiveDataTypeDefinitions.h"
#include "preProcessorDirectives.h"

using namespace std;

/* ------------------------------------------------------------------------- */
/* STATIC QUEUE/STACK CLASS */
/*
Objective:
	* Provide a Static Queue/Stack Implementation:
Rationale:
    * Since dynamic memory allocation can be expensive, we want to provide Queue
	  implementation that initializes data only once.
Assumption:
	* The maximum size of the number of elements that can ever be in the Queue is
	   know apriori.
	* Supports elements only of type Integer (regular or long)
	* Will add one extra element to the vector ( maxSize + 1 ) to wrap around
Functions Provided:
	  - Default Constructor : O(C)
	  - Constructor with given size : O(N)  (will have N+1 capacity, for wrap around)
	  - push_back(i) : O(C)
	  - front() : O(C)
	  - pop_front() : O(C) : !!! Modified, will return the element !!!
	  - empty() : O(C)
	  - clear() : O(C) //O(N): in regular case
	  - back()  : O(C)
	  - pop_back() : O(C) : !!! Modified, will return the element !!!
	  - size(): O(C)
*/
class staticQueue
{
	private:
		vector<MilanLongInt> squeue;
		MilanLongInt squeueHead;
		MilanLongInt squeueTail;
		MilanLongInt NumNodes;

		//Prevent Assignment and Pass by Value:
		staticQueue(const staticQueue& src);
		staticQueue& operator=(const staticQueue& rhs);

	public:
		//Constructors and Destructors
		staticQueue() { squeueHead = 0; squeueTail = 0; NumNodes = 0; }  //Default Constructor
		staticQueue(MilanLongInt maxSize) //MaximumSize
		{
			squeueHead = 0; //Head of the static Stack
			squeueTail = 0; //Tail of the Statuc Stack
			NumNodes = maxSize;
			try
			{
				squeue.reserve(NumNodes+1); //The number of nodes plus one to swap around
			}
			catch ( length_error )
			{
				cerr<<"Within Function: staticQueue(MilanLongInt maxSize) \n";
				cerr<<"Error: Not enough memory to allocate for Queue \n";
				exit(1);
			}
			squeue.resize( NumNodes+1, -1 ); //Initialize the stack with -1
		}
		~staticQueue() {};  //The destructor

		//Access:
		MilanLongInt front()  {	return squeue[squeueHead]; } //Non destructive
		MilanLongInt back()
		{
			if ( squeueTail == 0 ) //make it wrap around
				return squeue[NumNodes];
			else
				return squeue[squeueTail-1];
		}
		MilanLongInt getHead() { return squeueHead; }
		MilanLongInt getTail() { return squeueTail; }

		//Manipulation:
		void push_back(MilanLongInt newElement)
		{
		  //Q.push_back(i);
			squeue[squeueTail] = newElement;
			squeueTail = (squeueTail+1)%(NumNodes+1);
		}
		MilanLongInt pop_front() // !!! Modified, will return the element !!!
		{
		   //Q.pop_front();
			MilanLongInt U = squeue[squeueHead];
			squeueHead = (squeueHead+1)%(NumNodes+1);
			return U;
		}
		MilanLongInt pop_back() //!!! Modified, will return the element !!!
		{
			//S.pop_back();
			if ( squeueTail == 0 ) //make it wrap around
				squeueTail = NumNodes;
			else
				squeueTail = (squeueTail-1); //Remove the last element
			return squeue[squeueTail]; //Needs to be here. Because, the tail always points to the
			 //counter after the last existing element.
		}
		void clear()
		{
			//Q.clear(); //Empty the Queue
			squeueHead = 0; //Head of the static Queue
			squeueTail = 0; //Tail of the Statuc Queue
		}

		//Query:
		MilanBool empty()
		{
			//Q.empty();
			if ( squeueHead == squeueTail )
				return true;
			else
				return false;
		} //end of empty()
		MilanLongInt size()
		{
			//Q.size();
			MilanLongInt size = 0;
			if ( squeueHead == squeueTail )
				return size;
			else
				if ( squeueHead < squeueTail )
					return ( squeueTail - squeueHead );
				else
					return ( NumNodes + 1 - squeueHead + squeueTail );
		} //End of size()
		void display()
		{
			//Q.display();
			MilanLongInt i=0;
			cout<<"Queue: "<<endl;
			if ( squeueHead == squeueTail )
				cout<<"Empty"<<endl;
			else
				if ( squeueHead < squeueTail )
				{
					for ( i=squeueHead; i<squeueTail; i++ )
						cout<<squeue[i]<<", ";
					cout<<endl;
				}
				else
				{
					for ( i=squeueHead; i<NumNodes; i++)
						cout<<squeue[i]<<", ";
					for ( i=0; i<squeueTail; i++)
						cout<<squeue[i]<<", ";
					cout<<endl;
				}
		} //End of display()

};

#endif
