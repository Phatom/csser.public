/*
 * SfluidWindingObservable.cpp
 *
 *  Created on: May 30, 2014
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

SfluidWindingObservable::SfluidWindingObservable(string s, core::StateVariable* lstate,
core::Hamiltonian* lham, FILE* llog) : Measurable(s,lstate,lham,llog)
{
	Zcount = 0;
	wind2[0] = wind2[1] = wind2[2] = 0;
}

SfluidWindingObservable::~SfluidWindingObservable()
{

}

void inline SfluidWindingObservable::setup(string s, core::StateVariable* lstate, core::Hamiltonian* lham, FILE* llog)
{
	this->fileName = s;
	this->state = lstate;
	this->hamiltonian = lham;
	this->log = llog;

	Zcount = 0;
	wind2[0] = wind2[1] = wind2[2] = 0;
}

void inline SfluidWindingObservable::display()
{
	fprintf(log,"%15d\t %15d\t %15d",wind2[0],wind2[1],wind2[2]);
}

int inline SfluidWindingObservable::read()
{
	return 0;
}

int inline SfluidWindingObservable::write()
{
	ofstream rif(fileName.c_str(),ios::app);
	rif.precision(FIELDPRECISION);
	rif.width(FIELDWIDTH);
	rif.setf(FIELDFORMAT);

	double iZ = 1.0/Zcount;
	double factor = 0.5/this->state->beta*this->hamiltonian->nspd*this->hamiltonian->nspd;
	double w2x = wind2[0]*iZ*factor;
	double w2y = wind2[1]*iZ*factor;
	double w2z = wind2[2]*iZ*factor;
	double w2 = (w2x+w2y+w2z)/3;

	rif<<w2x<<" "<<w2y<<" "<<w2z<<" "<<w2<<endl;
	rif.close();

	wind2[0] = wind2[1] = wind2[2] = 0.0;
	Zcount = 0;

	return 0;
}

int inline SfluidWindingObservable::readViaIndex(int idx)
{

}

int inline SfluidWindingObservable::writeViaIndex(int idx)
{
	//Ignore idx
	write();
}

void inline SfluidWindingObservable::measure()
{
	int lwind[3] = {0,0,0};

	for(int i=0;i<state->offdiagOperatorLocs.size();i++)
	{
		unsigned int opstate = state->opString[state->offdiagOperatorLocs[i]];
		int posloc = opstate%hamiltonian->BT; //0 means diagonal 1,2 are off-diagonal entries
		unsigned int bond = opstate/hamiltonian->BT;
		int tt = (posloc == 1) ? 1 : -1;
		char bstate = hamiltonian->bondMap[bond*hamiltonian->bondMapOffset + 2];
		if(bstate!=0) //This is a edge state
			lwind[bstate-1] += tt;
	}

	Zcount++;
	this->wind2[0] += lwind[0]*lwind[0];
	this->wind2[1] += lwind[1]*lwind[1];
	this->wind2[2] += lwind[2]*lwind[2];
}

void inline SfluidWindingObservable::clear()
{
	wind2[0] = wind2[1] = wind2[2] = 0.0;
	Zcount = 0;
}

int inline SfluidWindingObservable::checkPointWrite()
{
	ofstream wif(chkFileName);
	if(wif)
		wif<<Zcount<<" "<<wind2[0]<<" "<<wind2[1]<<" "<<wind2[2]<<endl;
	wif.close();
	return 0;
}

int inline SfluidWindingObservable::checkPointRead()
{
	ifstream rif(chkFileName);
	if(rif)
		rif>>Zcount>>wind2[0]>>wind2[1]>>wind2[2];
	rif.close();
	return 0;
}


