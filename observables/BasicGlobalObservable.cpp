/*
 * BasicObservable.cpp
 *
 *  Created on: Apr 22, 2014
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

using namespace std;
using namespace __gnu_cxx;
using namespace measures;

BasicGlobalObservable::BasicGlobalObservable(string s, core::StateVariable* lstate,
core::Hamiltonian* lham, FILE* llog) : Measurable(s,lstate,lham,llog)
{
	Cbond = hamiltonian->C*hamiltonian->numBonds;

	staticO.operatorCount = 0;
	staticO.density = 0;
	staticO.energy = 0.0;
}

BasicGlobalObservable::~BasicGlobalObservable()
{
	this->collection.clear();
}

void inline BasicGlobalObservable::setup(string s, core::StateVariable* lstate, core::Hamiltonian* lham, FILE* llog)
{
	this->fileName = s;
	this->state = lstate;
	this->hamiltonian = lham;
	this->log = llog;

	Cbond = hamiltonian->C*hamiltonian->numBonds;

	staticO.operatorCount = 0;
	staticO.density = 0;
	staticO.energy = 0.0;
}

void inline BasicGlobalObservable::display()
{
	fprintf(log,"%15d\t %8d\t %10.6le ",staticO.operatorCount,staticO.density,staticO.energy);
}

int inline BasicGlobalObservable::read()
{
	unsigned int op, op2, rho, rho2;
	double e,e2;

	ifstream rif(fileName.c_str());
	while(rif)
	{
		rif>>op>>rho>>e;
		if(!rif)
			break;

		collection.push_back(staticObservable(op,rho,e));
	}
	rif.close();

	return 0;
}

int inline BasicGlobalObservable::write()
{
	ofstream rif(fileName.c_str(),ios::app);
	rif.precision(FIELDPRECISION);
	rif.width(FIELDWIDTH);
	rif.setf(FIELDFORMAT);

	for(vector<staticObservable>::iterator it=collection.begin();it!=collection.end();++it)
		rif<<it->operatorCount<<" "<<it->density<<" "<<it->energy<<endl;
	rif.close();

	return 0;
}

int inline BasicGlobalObservable::readViaIndex(int idx)
{
	stringstream ss;
	ss << idx;

	unsigned int op, op2, rho, rho2;
	double e,e2;

	ifstream rif(fileName.c_str() + ss.str());
	while(rif)
	{
		rif>>op>>rho>>e;
		if(!rif)
			break;

		collection.push_back(staticObservable(op,rho,e));
	}
	rif.close();

	return 0;
}

int inline BasicGlobalObservable::writeViaIndex(int idx)
{
	stringstream ss;
	ss << idx;

	ofstream rif(fileName.c_str() + ss.str(),ios::app);
	rif.precision(FIELDPRECISION);
	rif.width(FIELDWIDTH);
	rif.setf(FIELDFORMAT);

	for(vector<staticObservable>::iterator it=collection.begin();it!=collection.end();++it)
		rif<<it->operatorCount<<" "<<it->density<<" "<<it->energy<<endl;
	rif.close();

	clear();
	return 0;
}

void inline BasicGlobalObservable::measure()
{
	staticO.operatorCount = state->nOp;
	staticO.density = 0;
	for(int i=0;i<state->statesize;i++)
		staticO.density+=state->state[i];
	staticO.energy = -state->nOp*state->ibeta + Cbond + hamiltonian->mu*staticO.density;

	collection.push_back(staticObservable(staticO.operatorCount,staticO.density,staticO.energy));
}

void inline BasicGlobalObservable::clear()
{
	staticO.density = 0;
	staticO.energy = 0.0;
	staticO.operatorCount = 0;
	collection.clear();
}

int inline BasicGlobalObservable::checkPointWrite()
{
	ofstream rif(chkFileName.c_str(),ios::app);
	for(vector<staticObservable>::iterator it=collection.begin();it!=collection.end();++it)
		rif<<it->operatorCount<<" "<<it->density<<std::setprecision(9)<<it->energy<<endl;
	rif.close();

	return 0;
}

int inline BasicGlobalObservable::checkPointRead()
{
	unsigned int op, op2, rho, rho2;
	double e,e2;

	ifstream rif(chkFileName.c_str());
	while(rif)
	{
		rif>>op>>rho>>e;
		if(!rif)
			break;

		collection.push_back(staticObservable(op,rho,e));
	}
	rif.close();

	return 0;
}
