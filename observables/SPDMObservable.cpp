/*
 * SPDMDynaMeasurable.cpp
 *
 *  Created on: Apr 29, 2014
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

SPDMObservable::SPDMObservable(string s, core::StateVariable* lstate,
core::Hamiltonian* lham, FILE* llog, map2d& _spdm) : Measurable(s,lstate,lham,llog), spdm(_spdm)
{
	Zcount = 0;
}

SPDMObservable::~SPDMObservable()
{
	Zcount = 0;
	clear();
}

void SPDMObservable::measure()
{
	Zcount++;
}

void SPDMObservable::clear()
{
	Zcount = 0;
	map2d::iterator it = spdm.begin();
	for(;it!=spdm.end();++it)
		it->second.clear();
	spdm.clear();
}

void SPDMObservable::display()
{

}

int SPDMObservable::write()
{

}

int SPDMObservable::read()
{

}

int SPDMObservable::writeViaIndex(int idx)
{
	double iZ = 1.0/Zcount/state->lpContrIter;
	Zcount = 0;

	stringstream ss;
	ss<<idx;
	string tfname = this->fileName + ss.str();

	ofstream wfp(tfname);
	wfp.precision(FIELDPRECISION);
	wfp.width(FIELDWIDTH);
	wfp.setf(FIELDFORMAT);

	map2d::iterator it = spdm.begin();
	for(;it!=spdm.end();++it)
	{
		for(Row::iterator it2 = it->second.begin();it2!=it->second.end();++it2)
		{
			double mval = iZ*it2->second;
			wfp<<it->first<<" "<<it2->first<<" "<<mval<<endl;
		}

		it->second.clear();
	}
	wfp.close();
	spdm.clear();
}

int SPDMObservable::readViaIndex(int idx)
{
	stringstream ss;
	ss<<idx;
	string tfname = this->fileName + ss.str();

	int i,j;
	double val;

	ifstream rfp(tfname);
	while(rfp)
	{
		rfp>>i>>j>>val;

		spdm[i][j] = val;
	}
	rfp.close();
}

int SPDMObservable::checkPointWrite()
{
	ofstream wfp(chkFileName);
	wfp.precision(FIELDPRECISION);
	wfp.width(FIELDWIDTH);
	wfp.setf(FIELDFORMAT);

	map2d::iterator it = spdm.begin();
	for(;it!=spdm.end();++it)
		for(Row::iterator it2 = it->second.begin();it2!=it->second.end();++it2)
			wfp<<it->first<<" "<<it2->first<<" "<<it2->second<<endl;

	wfp.close();
}

int SPDMObservable::checkPointRead()
{
	int i,j;
	double val;

	ifstream rfp(chkFileName);
	while(rfp)
	{
		rfp>>i>>j>>val;
		spdm[i][j] = val;
	}
	rfp.close();
}
/* namespace measures */
