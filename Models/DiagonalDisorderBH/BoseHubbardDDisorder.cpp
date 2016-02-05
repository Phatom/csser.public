/*
 * BoseHubbardDDisorder.cpp
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
#include "BoseHubbardDDisorder.h"

using namespace std;
using namespace __gnu_cxx;
using namespace core;

BoseHubbardDDisorder::BoseHubbardDDisorder(int _procId, string _path, string fileName, FILE* llog) {
	// TODO Auto-generated constructor stub
	this->processId = _procId;
	this->path = _path;

	this->latticeParamFileName = fileName;
	this->BT = 3;
	this->bondMapOffset = 3;
	this->symmSetOffset = 3;
	log=llog;

	rgenref = gsl_rng_alloc(gsl_rng_mt19937);
	disorderFileName = "defaultDisorderFile.dat";
}

BoseHubbardDDisorder::~BoseHubbardDDisorder()
{
	delete[] this->bondMap;
	delete[] this->bondToSSet;
	delete[] this->E;
}

//Recall 0 = a+, 1 = a
double inline BoseHubbardDDisorder::operateWith(char op,char state)
{
	return sqrt(state + 1-op);
}

double inline BoseHubbardDDisorder::getDiagonalTerm(int si,int sj,int statei, int statej)
{
	return (this->C + multA*(-(E[si]*statei + E[sj]*statej) - 0.5*(U*statei*(statei - 1) + U*statej*(statej - 1))));
}

double inline BoseHubbardDDisorder::getOffdiagonalTerm(unsigned int bond, char op, const char* state)
{
	char b0 = state[2]-state[0];
	char b1 = state[3]-state[1];

	b0 = (b0<0) ? 0:1;
	b1 = (b1<0) ? 0:1;

	return t*sqrt((state[0]+b0)*(state[1]+b1));
}

int BoseHubbardDDisorder::loadHamiltonianParams()
{
	string temp;
	float trapFreq;

	ifstream f(this->latticeParamFileName.c_str());
	if(f)
	{
		f>>t>>U;
		f>>mu;
		f>>Delta;
		f>>trapFreq;
		this->Omega = 0.5*87*(1.66053886e-27/2.31e-30)*406*406*1.0e-18*4*(9.869604401)*trapFreq*trapFreq; //0.5*m*a^2*(2*pi)^2/Er*freq^2

		f>>this->dimension;
		f>>this->nspd;
		f>>this->LAM;
		f>>this->EPS;

		f>>temp;
		this->bc = (temp.compare("pbc")==0) ? PBC: HW;

		f>>temp;
		f>>disorderSeed;
		if(temp.compare("True")==0) //Randomize disorder
			disorderSeed = (int) time(NULL)%(65535);

		f>>temp;
		if(temp.compare("True")==0) //Apply process id
			disorderSeed += processId;

		f>>disorderFileName; disorderFileName = path + disorderFileName;
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

	//Create diagonal disorder and trap. Also subsume chemical potential into the diagonal term
	this->E = new double[this->numSites];

	ifstream disf(disorderFileName);
	if(disf)
	{
		fprintf(log,"Disorder loaded from file %s\n",disorderFileName.c_str());

		int i;
		for(i=0;i<numSites;i++)
		{
			disf>>E[i];
			if(!disf)
				break;
		}
		disf.close();

		if(i!=numSites)
		{
			fprintf(log,"Disorder file data cuts-off abruptly. Size should be %d, but current size is %d\n",numSites,i);
			return -1;
		}
	}
	else
	{
		fprintf(log,"Generating new disorder\n");
		gsl_rng_set(rgenref,disorderSeed);
		for(int i=0;i<numSites;i++)
			E[i] = this->Delta*(gsl_rng_uniform(rgenref)-0.5)*2.0;

		ofstream disf(disorderFileName);
		for(int i=0;i<numSites;i++)
			disf<<E[i]<<endl;
		disf.close();
	}

	int X,Y,Z;
	X=Y=Z=this->nspd/2; //it's great if it's odd otherwise it wont' be nice and symmetric

	if(this->dimension == 1)
	{
		for(int x=0;x<nspd;x++)
		{
			double xpow = pow(x-X,2);
			E[x] += Omega*xpow - mu;
		}
	}
	else if(this->dimension == 2)
	{
		for(int x=0;x<nspd;x++)
			for(int y=0;y<nspd;y++)
			{
				double xpow = pow(x-X,2);
				double ypow = pow(y-Y,2);
				E[nspd*y + x] += Omega*(xpow+ypow) - mu;
			}
	}
	else if(this->dimension == 3)
	{
		int ns2 = nspd*nspd;
		for(int x=0;x<nspd;x++)
			for(int y=0;y<nspd;y++)
				for(int z=0;z<nspd;z++)
				{
					double xpow = pow(x-X,2);
					double ypow = pow(y-Y,2);
					double zpow = pow(z-Z,2);
					E[ns2*z + nspd*y + x] += Omega*(xpow+ypow+zpow) - mu;
				}
	}
	else
	{
		fprintf(log,"Dimension: %d not supported\n",dimension);
		return -1;
	}

}

int BoseHubbardDDisorder::createSymmetricSet()
{
	//Use the bondMap directly
	//We are writing this is in the form of (si,sj)
	this->symmSetOffset = 3;
	this->symmetricSet = this->bondMap+3; //Skip 0 since it is reserved for unit operators
	this->symmetricSetSize = this->numBonds;

	//Do the bond to symmetric set assignment
	for(int i=1;i<=numBonds;i++)
		bondToSSet[i] = i-1;
}

void BoseHubbardDDisorder::calculateConstant()
{
	C=0.0;
	double tv;

	//Calculate constants
	for(long b = 1;b<=numBonds;b++)
	{
		int si = bondMap[this->bondMapOffset*b];
		int sj = bondMap[this->bondMapOffset*b + 1];

		for(char ni=0;ni<(LAM+1);ni++)
			for(char nj=0;nj<(LAM+1);nj++)
			{
				tv = multA*( -(ni*E[si]+nj*E[sj]) - 0.5*(U*ni*(ni-1) + U*nj*(nj-1)));

				if(tv<0 && tv<C)
					C = tv;
			}
	}

	C = -C;
	C += EPS;
}

void BoseHubbardDDisorder::display()
{
	fprintf(log,"==========================================================================\n");
	fprintf(log,"Hamiltonia Parameters\n");
	fprintf(log,"==========================================================================\n");
	fprintf(log,"Boundary Conditions: %s\n",((this->bc == PBC) ? "PBC": "HW"));
	fprintf(log,"t: %9.3e\n",this->t);
	fprintf(log,"U: %9.3e\n",this->U);
	fprintf(log,"mu: %16.10e\n",this->mu);
	fprintf(log,"Disorder Bound: %16.10e\n",this->Delta);
	fprintf(log,"Trap Curvature: %16.10e\n",this->Omega);
	fprintf(log,"Disorder Seed: %d\n",this->disorderSeed);

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
