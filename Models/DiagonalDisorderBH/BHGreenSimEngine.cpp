/*
 * BHGreenSimEngine.cpp
 *
 *  Created on: May 5, 2014
 *      Author: ushnish
 *
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
#include "BHGreenSimEngine.h"

namespace measures {

int inline BHGreenSimEngine::sweep(int& duper, long long& luper)
{
	do{

		duper = this->diagonalUpdate();

		//Expand if possible
		if(duper==-1)
		{
			int ecode = this->stateVariable->reSize(1.25*this->stateVariable->stringsize);
			if(ecode==EXITCODE)
				return EXITCODE;
			this->stateVariable->reset();
		}
	}while(duper==-1);

	this->maxLoopLength = 100*this->stateVariable->nOp;

	int ecode = this->loopUpdate(luper);

	if(ecode == EXITCODE)
		return EXITCODE;
	else
		return 0;
}

} /* namespace measures */
