/*
 * BasicStaticObservable.cpp
 *
 *  Created on: Aug 8, 2014
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

namespace measures {

BasicStaticObservable::BasicStaticObservable(string s, core::StateVariable* lstate,
core::Hamiltonian* lham, FILE* llog) : Measurable(s,lstate,lham,llog)
{
	Cbond = hamiltonian->C*hamiltonian->numBonds;
	operatorCount = 0;
	density = 0;
	energy = 0.0;
}

BasicStaticObservable::~BasicStaticObservable()
{

}

void inline BasicStaticObservable::setup(string s, core::StateVariable* lstate, core::Hamiltonian* lham, FILE* llog)
{
	this->fileName = s;
	this->state = lstate;
	this->hamiltonian = lham;
	this->log = llog;

	Cbond = hamiltonian->C*hamiltonian->numBonds;

	operatorCount = 0;
	density = 0;
	energy = 0.0;
}

void inline BasicStaticObservable::display()
{
	fprintf(log,"%15d\t %8d\t %10.6le ",operatorCount,density,energy);
}

int inline BasicStaticObservable::read()
{
	return 0;
}

int inline BasicStaticObservable::write()
{
	return 0;
}

int inline BasicStaticObservable::readViaIndex(int idx)
{
	return 0;
}

int inline BasicStaticObservable::writeViaIndex(int idx)
{
	return 0;
}

void inline BasicStaticObservable::measure()
{
	operatorCount = state->nOp;
	density = 0;
	for(int i=0;i<state->statesize;i++)
		density+=state->state[i];
	energy = -state->nOp*state->ibeta + Cbond + hamiltonian->mu*density;
}

void inline BasicStaticObservable::clear()
{
	density = 0;
	energy = 0.0;
	operatorCount = 0;
}

int inline BasicStaticObservable::checkPointWrite()
{
	return 0;
}

int inline BasicStaticObservable::checkPointRead()
{
}

} /* namespace measures */
