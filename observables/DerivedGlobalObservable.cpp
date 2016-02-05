/*
 * DerivedGlobalObservable.cpp
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

DerivedGlobalObservable::DerivedGlobalObservable(string s, core::StateVariable* lstate, core::Hamiltonian* lham, FILE* llog) : Measurable(s,lstate,lham,llog)
{
	sample = 0;
	opCountavg = 0;
	opCount2avg = 0;
	densityavg = 0;
	density2avg = 0;
}

DerivedGlobalObservable::~DerivedGlobalObservable()
{

}

void inline DerivedGlobalObservable::display()
{
	//Can be defined later if needed
}

int inline DerivedGlobalObservable::read()
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

int inline DerivedGlobalObservable::write()
{
	double is = 1.0/sample;
	ofstream rif(fileName.c_str(),ios::app);
	rif.precision(FIELDPRECISION);
	rif.width(FIELDWIDTH);
	rif.setf(FIELDFORMAT);

	double rho = densityavg*is;
	double rho2 = density2avg*is;
	double kappa = (rho2-rho*rho)*state->beta;

	double opavg = opCountavg*is;
	double op2avg = opCount2avg*is;
	double heatcap = op2avg - opavg*opavg - opavg;

	rif<<kappa<<" "<<heatcap<<endl;
	rif.close();

	sample = 0;
}

int inline DerivedGlobalObservable::readViaIndex(int idx)
{

}

int inline DerivedGlobalObservable::writeViaIndex(int bin)
{
	double is = 1.0/sample;
	ofstream rif(fileName.c_str(),ios::app);
	rif.precision(FIELDPRECISION);
	rif.width(FIELDWIDTH);
	rif.setf(FIELDFORMAT);

	double rho = densityavg*is;
	double rho2 = density2avg*is;
	double kappa = (rho2-rho*rho)*state->beta;

	double opavg = opCountavg*is;
	double op2avg = opCount2avg*is;
	double heatcap = op2avg - opavg*opavg - opavg;

	rif<<kappa<<" "<<heatcap<<endl;
	rif.close();

	sample = 0;
	opCountavg = 0;
	opCount2avg = 0;
	densityavg = 0;
	density2avg = 0;
}

void inline DerivedGlobalObservable::measure()
{
	sample++;
	opCountavg += state->nOp;
	opCount2avg += state->nOp*state->nOp;

	unsigned int totalrho = 0;
	for(int i=0;i<state->statesize;i++)
		totalrho += state->state[i];

	densityavg += totalrho;
	density2avg += totalrho*totalrho;

}

void inline DerivedGlobalObservable::clear()
{
	sample = 0;
	opCountavg = 0;
	opCount2avg = 0;
	densityavg = 0;
	density2avg = 0;
}

int inline DerivedGlobalObservable::checkPointWrite()
{
	ofstream rif(chkFileName.c_str());
	rif<<sample<<" "<<opCountavg<<" "<<opCount2avg<<" "<<densityavg<<" "<<density2avg<<endl;
	rif.close();
}

int inline DerivedGlobalObservable::checkPointRead()
{
	double lr,lr2;
	int p = 0;
	ifstream rif(chkFileName.c_str());
	rif>>sample>>opCountavg>>opCount2avg>>densityavg>>density2avg;
	rif.close();
}


