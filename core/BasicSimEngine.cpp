/*
 * BasicSimEngine.cpp
 *
 *  Created on: Apr 21, 2014
 *      Author: ushnish

  	Copyright (c) 2014 Ushnish Ray
	All rights reserved.

	This program and the accompanying materials are made available under explicit agreement
	between Ushnish Ray and the end user. You may not redistribute the code and
	accompanying material to anyone.

	On the event that the software is used to generate data that is used implicitly or explicitly
	for research purposes, proper acknowledgment must provided  in the citations section of
	publications.

	This software is cannot be used for commercial purposes in any way whatsoever.
 */

#include "mainincs.h"
#include "core.h"
#include "BasicSimEngine.h"

using namespace std;
using namespace __gnu_cxx;
using namespace core;

int BasicSimEngine::diagonalUpdate()
{
	//Do a basic calculation so that it is usable in loop
	multB = stateVariable->beta*hamiltonian->numBonds;
	imultB = 1.0/multB;

	int flag;

	memcpy(state_t,stateVariable->state,hamiltonian->numSites*sizeof(char));

	for(unsigned long j=0;j<stateVariable->opCutOff;j++)
	{

		char optype = stateVariable->opString[j]%hamiltonian->BT;

		if(stateVariable->opString[j] == 0) //Unit operator
		{

			unsigned int b = gsl_rng_uniform_int(rgenref,hamiltonian->numBonds)+1;

			unsigned int boffset = hamiltonian->bondMapOffset*b;
			unsigned int b_i = hamiltonian->bondMap[boffset];
			unsigned int b_j = hamiltonian->bondMap[boffset+1];
			char si = state_t[b_i];
			char sj = state_t[b_j];
			double dgProb = hamiltonian->getDiagonalTerm(b_i,b_j,si,sj);

			/*
			 (
			 global.C + multA*(global.mu*(si + sj)
			 - (global.EE[b_i]*si + global.EE[b_j]*sj)
			 - 0.5*(global.U[b_i]*si*(si - 1)
			 + global.U[b_j]*sj*(sj - 1))) - global.V*si*sj
			 );
			*/
			double trPb = multB*dgProb/(stateVariable->opCutOff - stateVariable->nOp);
			double rnd = gsl_rng_uniform(rgenref);
			if(rnd<=trPb)
			{

				stateVariable->opString[j] = hamiltonian->BT*b;
				stateVariable->nOp = stateVariable->nOp + 1;

				if(stateVariable->nOp >= stateVariable->opCutOff)
				{
					flag = 1; //Operator Cut Off changed!
					stateVariable->opCutOff = stateVariable->opCutOff*1.3;
					while(stateVariable->opCutOff <= stateVariable->nOp)
						stateVariable->opCutOff = stateVariable->opCutOff + 1;

#ifdef BOUNDCHECK
					if(stateVariable->opCutOff >= stateVariable->stringsize)
					{
						return -1;
					}
#endif
				}
			}
		}
		else if(optype == 0) //Even bond number so diagonal term
		{

			unsigned int b = stateVariable->opString[j]/hamiltonian->BT;
			unsigned int boffset = b*hamiltonian->bondMapOffset;

			unsigned int b_i = hamiltonian->bondMap[boffset];
			unsigned int b_j = hamiltonian->bondMap[boffset+1];
			char si = state_t[b_i];
			char sj = state_t[b_j];
			double dgProb = hamiltonian->getDiagonalTerm(b_i,b_j,si,sj);;

			/*
			(global.C + multA*(global.mu*(si + sj)
			- (global.EE[b_i]*si + global.EE[b_j]*sj)
			- 0.5*(global.U[b_i]*si*(si - 1)
			+ global.U[b_j]*sj*(sj - 1))) - global.V*si*sj);
			 */

			double trPb = 1.0*(stateVariable->opCutOff - stateVariable->nOp + 1)*imultB/dgProb;
			double rnd = gsl_rng_uniform(rgenref);
			if(rnd<=trPb)
			{
				stateVariable->opString[j] = 0;
				stateVariable->nOp = stateVariable->nOp - 1;
			}
		}
		else if(optype != 0) //Odd bond number so off-diagonal term
		{
			unsigned int b = stateVariable->opString[j]/hamiltonian->BT;
			unsigned int boffset = hamiltonian->bondMapOffset*b;
			//Propagate  opType = 1 => a(^+)a, opType = 2 => aa(^+).
			char uvi = (optype == 1) ? 1: -1;
			state_t[hamiltonian->bondMap[boffset]] += uvi;
			state_t[hamiltonian->bondMap[boffset+1]] -= uvi;
		}
	}
}


