/*
 * EquibRunMode.cpp
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
#include "EquibRunMode.h"

using namespace std;
using namespace __gnu_cxx;
using namespace runmode;
using namespace measures;

EquibRunMode::~EquibRunMode()
{
	delete eqbBase;
}

EquibRunMode::EquibRunMode(core::Hamiltonian* h,core::TransitionMap* t,core::StateVariable* s,core::SimulationEngine* seng,
		EqbParameters& ep, Flags* f, FILE* llog):RunMode(h,t,s,seng)
{
	eparams = ep;
	this->flag = f;
	log = llog;
	eqbBase = new BasicGlobalObservable("eqbState",s,h,llog);

	currState.beta = ep.startBeta;
	currState.sweep = 0;
}

int EquibRunMode::loadEqbState()
{
	ifstream rif(eparams.eqbStateFile);
	if(rif)
	{
		rif>>currState.beta>>currState.sweep;
		rif.close();
	}
	else
		return FILENOTFOUND;

	int ret = this->runState->readState();
	return ret;
}

void EquibRunMode::saveEqbState()
{
	ofstream oif(eparams.eqbStateFile);
	oif<<currState.beta<<" "<<currState.sweep;
	oif.close();

	this->runState->writeState();
}

void EquibRunMode::saveCheckPoint()
{
	string chkFileName = eparams.eqbStateFile + ".chk";
	ofstream oif(chkFileName);
	oif<<currState.beta<<" "<<currState.sweep;
	oif.close();

	string scfilename = this->runState->stateFileName + ".chk";
	this->runState->writeState(scfilename);
}

void EquibRunMode::saveTemperatureState(int foffset)
{
	stringstream sfo;
	sfo<<foffset;
	string tfile = this->runState->stateFileName + "_" + sfo.str();
	this->runState->writeState(tfile);
}
