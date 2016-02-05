/*
 * checkReadState.cpp
 *
 *  Created on: Jul 20, 2014
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

#include "sse.h"
#include "StateVariableBH.h"
#include "BoseHubbardDDisorder.h"
#include "EqbSimEngine.h"
#include "BHspdmSimEngine.h"
#include "BHGreenSimEngine.h"

using namespace std;
using namespace core;
using namespace measures;
using namespace runmode;

int main(int argc, char* argv[])
{
	if(argc<5)
	{
		cout<<"Format is: <State Parameter File> <State File> <output> <Output File>\n";
		return -1;
	}

	StateVariableBH* svar;
	string stateParamFile = string(argv[1]);
	string stateOutputFile = string(argv[2]);

	FILE* stdo = fopen(argv[3],"w");
	try{
		svar = new StateVariableBH(stdo,stateParamFile,stateOutputFile);
	}
	catch(int exception)
	{
		if(exception == ConstructionException)
		{
			fprintf(stdo,"String Size is larger than HARDCUTOFF. Re-compile after increasing size.\n");
			fflush(stdo);
			return EXITCODE;
		}
	}


	svar->readState();
	svar->display(stdo);

	FILE* wf = fopen(argv[4],"w");
	for(int i=0;i<svar->statesize;i++)
	{
		fprintf(wf,"%d\n",svar->state[i]);
	}
	fclose(wf);

	fclose(stdo);
	delete svar;
	return 0;
}



