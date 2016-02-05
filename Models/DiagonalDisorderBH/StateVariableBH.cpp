/*
 * StateVariableBH.cpp
 *
 *  Created on: Apr 24, 2014
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
#include "StateVariableBH.h"

using namespace std;

StateVariableBH::StateVariableBH(FILE* _stdo, string startStateFile, string saveLocation)
{
	this->stdo = _stdo;

	ifstream rif(startStateFile.c_str());
	if(rif)
	{
		rif>>this->LAM;
		rif>>this->statesize;
		rif>>this->opCutOff;
		rif>>this->stringsize;
		rif>>this->lpContrIter;
	}
	rif.close();

	//Check on string size
	if(this->stringsize>HARDCUTOFF)
		throw ConstructionException;

	//Intentionally keeping this separate
	this->stateFileName = saveLocation;
	this->allocate();
}

int StateVariableBH::initializeState()
{
	//Random initialization
	for(int i=0;i<this->statesize;i++)
		this->state[i] = gsl_rng_uniform_int(rgenref,LAM+1);
}
