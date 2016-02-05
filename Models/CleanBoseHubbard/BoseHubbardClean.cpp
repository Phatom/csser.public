/*
 * BoseHubbardClean.cpp
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
#include "BoseHubbardClean.h"

using namespace std;
using namespace __gnu_cxx;
using namespace core;

BoseHubbardClean::BoseHubbardClean(string fileName, FILE* llog) {
	// TODO Auto-generated constructor stub
	this->latticeParamFileName = fileName;
	this->BT = 3;
	this->bondMapOffset = 3;
	this->symmSetOffset = 2;
	log=llog;
}

BoseHubbardClean::~BoseHubbardClean()
{
	delete[] this->symmetricSet;
	delete[] this->bondMap;
	delete[] this->bondToSSet;
}

//Recall 0 = a+, 1 = a
double inline BoseHubbardClean::operateWith(char op,char state)
{
	return sqrt(state + 1-op);
}

double inline BoseHubbardClean::getDiagonalTerm(int si,int sj,int statei, int statej)
{
	return (this->C + multA*(mu*(statei + statej) - 0.5*(U*statei*(statei - 1) + U*statej*(statej - 1))));
}

double inline BoseHubbardClean::getOffdiagonalTerm(unsigned int bond, char op, const char* state)
{
	char b0 = state[2]-state[0];
	char b1 = state[3]-state[1];

	b0 = (b0<0) ? 0:1;
	b1 = (b1<0) ? 0:1;

	return t*sqrt((state[0]+b0)*(state[1]+b1));
}

int BoseHubbardClean::loadHamiltonianParams()
{
	string temp;

	ifstream f(this->latticeParamFileName.c_str());
	if(f)
	{
		f>>t>>U;
		f>>mu;
		f>>this->dimension;
		f>>this->nspd;
		f>>this->LAM;
		f>>this->EPS;

		f>>temp;
		this->bc = (temp.compare("pbc")==0) ? PBC: HW;
	}
	else
		return -1;

	f.close();

	this->numSites = pow(this->nspd,this->dimension);
	this->numBonds = numSites*this->dimension;
	if(this->bc==HW)
		this->numBonds -= this->dimension*pow(this->nspd,this->dimension-1);

	//Need one for spot other than number of bonds for unit operator
	bondMap = new unsigned int[(this->numBonds+1)*this->bondMapOffset];

	//Also need bond to symmetric set
	//Need to account for b = 0 which is the unit operator and doesn't belong to set
	this->bondToSSet = new unsigned int[this->numBonds+1];

	multA = 1.0/this->dimension*0.5;
}

int BoseHubbardClean::createSymmetricSet()
{
	//We are writing this is in the form of (si,sj)
	this->symmSetOffset = 2;
	this->symmetricSet = new unsigned int[2];
	this->symmetricSet[0] = 0;
	this->symmetricSet[1] = 1;
	this->symmetricSetSize = 1;

	//Do the bond to symmetric set assignment
	for(int i=1;i<=numBonds;i++)
		bondToSSet[i] = 0;
}

void BoseHubbardClean::calculateConstant()
{
	C=0.0;
	double tv;

	//Calculate constants
	for(long b = 1;b<=numBonds;b++)
		for(char ni=0;ni<(LAM+1);ni++)
			for(char nj=0;nj<(LAM+1);nj++)
			{
				//int si = bondMap[this->bondMapOffset*b];
				//int sj = bondMap[this->bondMapOffset*b + 1];

				tv = multA*( mu*(ni+nj) - 0.5*(U*ni*(ni-1) + U*nj*(nj-1)));

				if(tv<0 && tv<C)
					C = tv;
			}

	C = -C;
	C += EPS;
}

void BoseHubbardClean::display()
{
	fprintf(log,"==========================================================================\n");
	fprintf(log,"Hamiltonia Parameters\n");
	fprintf(log,"==========================================================================\n");
	fprintf(log,"Boundary Conditions: %s\n",((this->bc == PBC) ? "PBC": "HW"));
	fprintf(log,"t: %9.3e\n",this->t);
	fprintf(log,"U: %9.3e\n",this->U);
	fprintf(log,"mu: %16.10e\n",this->mu);
	fprintf(log,"LAM: %d\n",this->LAM);
	fprintf(log,"EPS: %9.3e\n",this->EPS);
	fprintf(log,"C: %10.6le\n",this->C);
	fprintf(log,"Lattice dimension: %u\n",this->dimension);
	fprintf(log,"Lattice linear dimension: %u\n",this->nspd);
	fprintf(log,"Number of Sites in Lattice: %u\n",this->numSites);
	fprintf(log,"Number of Bonds: %Ld\n",this->numBonds);
	fprintf(log,"\n\n");
	fflush(log);
}
