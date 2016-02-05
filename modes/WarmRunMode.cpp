/*
 * WarmRunMode.cpp
 *
 *  Created on: Apr 23, 2014
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
#include "observables.h"
#include "RunMode.h"
#include "WarmRunMode.h"

using namespace std;
using namespace __gnu_cxx;
using namespace runmode;
using namespace measures;

WarmRunMode::~WarmRunMode()
{
	delete eqbBase;
}

WarmRunMode::WarmRunMode(core::Hamiltonian* h,core::TransitionMap* t,core::StateVariable* s,core::SimulationEngine* seng,
		WarmParameters& ep, Flags* f, FILE* llog):RunMode(h,t,s,seng)
{
	params = ep;
	this->flag = f;
	log = llog;
	eqbBase = new BasicGlobalObservable("warmState",s,h,llog);

	currState.sweep = 0;
}

int WarmRunMode::loadState()
{
	bool pload = false;
	ifstream rif(params.warmStateFile);
	if(rif)
	{
		rif>>currState.sweep;
		rif.close();
	}
	else
		pload = true;

	int ret = this->runState->readState();

	if(pload && ret==FILENOTFOUND)
		return ret;
	else if(pload && ret==EXITCODE)
		return ret;
	else if(pload)
		return PARTIALFAIL;

	return ret;
}

void WarmRunMode::saveState()
{
	ofstream oif(params.warmStateFile);
	oif<<currState.sweep;
	oif.close();

	this->runState->writeState();
}

void WarmRunMode::saveCheckPoint()
{
	string chkFileName = params.warmStateFile + ".chk";
	ofstream oif(chkFileName);
	oif<<currState.sweep;
	oif.close();

	string scfilename = this->runState->stateFileName + ".chk";
	this->runState->writeState(scfilename);
}
