/*
 * BasicLocalObservable.cpp
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

using namespace std;
using namespace __gnu_cxx;
using namespace measures;

BasicLocalObservable::BasicLocalObservable(string s, core::StateVariable* lstate, core::Hamiltonian* lham, FILE* llog) : Measurable(s,lstate,lham,llog)
{
	sample = 0;
	localRho = new float[lstate->statesize];
	localRho2 = new float[lstate->statesize];
	for(int i=0;i<lstate->statesize;i++)
		localRho[i] = localRho2[i] = 0.0;
}

BasicLocalObservable::~BasicLocalObservable()
{
	delete[] localRho;
	delete[] localRho2;
}

void inline BasicLocalObservable::display()
{
	//Can be defined later if needed
}

int inline BasicLocalObservable::read()
{
	/*
	double lr,lr2;
	int p = 0;
	ifstream rif(fileName.c_str());
	rif>>sample;
	while(rif)
	{
		rif>>lr>>lr2;
		if(!rif)
			break;

		localRho[p] = lr;
		localRho2[p++] = lr2;
	}
	rif.close();
	*/
}

int inline BasicLocalObservable::readViaIndex(int idx)
{
	/*
	double lr,lr2;
	int p = 0;
	ifstream rif(fileName.c_str());
	rif>>sample;
	while(rif)
	{
		rif>>lr>>lr2;
		if(!rif)
			break;

		localRho[p] = lr;
		localRho2[p++] = lr2;
	}
	rif.close();
	*/
}

int inline BasicLocalObservable::write()
{
	//This time reset
	double is = 1.0/sample;
	ofstream rif(fileName);
	rif.precision(FIELDPRECISION);
	rif.width(FIELDWIDTH);
	rif.setf(FIELDFORMAT);

	for(int p=0;p<state->statesize;p++)
	{
		double rho = localRho[p]*is; localRho[p] = 0;
		double rho2 = localRho2[p]*is; localRho2[p] = 0;
		double kappa = (rho2-rho*rho)*state->beta;

		rif<<rho<<" "<<kappa<<endl;
	}
	rif.close();

	sample = 0;
	bin++;
}


int inline BasicLocalObservable::writeViaIndex(int idx)
{
	//This time reset
	stringstream ss;
	ss<<idx;
	double is = 1.0/sample;
	ofstream rif((fileName + "_" + ss.str()).c_str());
	rif.precision(FIELDPRECISION);
	rif.width(FIELDWIDTH);
	rif.setf(FIELDFORMAT);

	for(int p=0;p<state->statesize;p++)
	{
		double rho = localRho[p]*is; localRho[p] = 0;
		double rho2 = localRho2[p]*is; localRho2[p] = 0;
		double kappa = (rho2-rho*rho)*state->beta;

		rif<<rho<<" "<<kappa<<endl;
	}
	rif.close();

	sample = 0;
	bin++;
}


void inline BasicLocalObservable::measure()
{
	sample++;
	for(int i=0;i<state->statesize;i++)
	{
		localRho[i] += state->state[i];
		localRho2[i] += state->state[i]*state->state[i];
	}
}

void inline BasicLocalObservable::clear()
{
	bin = 0;
	sample = 0;
	for(int i=0;i<state->statesize;i++)
		localRho[i] = localRho2[i] = 0;
}

int inline BasicLocalObservable::checkPointWrite()
{
	ofstream rif(chkFileName.c_str());
	rif.precision(FIELDPRECISION);
	rif.width(FIELDWIDTH);
	rif.setf(FIELDFORMAT);

	rif<<bin<<endl;
	rif<<sample<<endl;
	for(int p=0;p<state->statesize;p++)
		rif<<localRho[p]<<" "<<localRho2[p]<<endl;
	rif.close();
}

int inline BasicLocalObservable::checkPointRead()
{
	double lr,lr2;
	int p = 0;
	ifstream rif(chkFileName.c_str());
	rif>>bin;
	rif>>sample;
	while(rif)
	{
		rif>>lr>>lr2;
		if(!rif)
			break;

		localRho[p] = lr;
		localRho2[p++] = lr2;
	}
	rif.close();
}


