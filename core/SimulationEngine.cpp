/*
 * SimulationEngine.cpp
 *
 *  Created on: Apr 21, 2014
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

using namespace std;
using namespace __gnu_cxx;
using namespace core;

SimulationEngine::SimulationEngine(Hamiltonian* h, TransitionMap* t, StateVariable* s, gsl_rng* _lrng, FILE* lstdo)
{
	hamiltonian = h;
	tmap = t;
	stateVariable = s;

	first = new long[hamiltonian->numSites];
	last = new long[hamiltonian->numSites];
	state_t = new char[hamiltonian->numSites];

	rgenref = _lrng;

	invNumBonds = 1.0/hamiltonian->numBonds;
	multB = stateVariable->beta*hamiltonian->numBonds;
	imultB = 1.0/multB;

	stdo = lstdo;
}

SimulationEngine::~SimulationEngine()
{
	delete[] first;
	delete[] last;
	delete[] state_t;
}

void SimulationEngine::setRandomGenerator(gsl_rng* _lgen)
{
	rgenref = _lgen;
}

void SimulationEngine::setLog(FILE* log)
{
	stdo = log;
}
