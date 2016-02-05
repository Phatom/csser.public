/*
 * stateVariable.cpp
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
#include "StateVariable.h"

using namespace std;
using namespace __gnu_cxx;
using namespace core;

StateVariable::StateVariable()
{
	this->allocated = false;
	this->beta = 0.0;
	this->lpContrIter = 0;
	this->nOp = 0;
	this->opCutOff = 0;
	this->statesize = 0;
	this->stringsize = 0;
	this->stateFileName = "";
}

StateVariable::StateVariable(int _statesize,int _strsize, long _opCutOff, FILE* _stdo, float _beta, int _lpc, string _stateFileName)
{
	stdo = _stdo;
	nOp = 0;
	beta = _beta; ibeta = 1.0/_beta;
	lpContrIter = _lpc;

	statesize = _statesize;
	stringsize = _strsize;
	opCutOff = _opCutOff;

	state = new char[statesize];
	opString = new bigi[stringsize];
	vtx = new bigi[stringsize];
	link = new bigi[4*stringsize];

	for(long i = 0; i<stringsize; i++)
	{
		opString[i] = 0;
		vtx[i] = 0;
	}

	for(long i=0; i<4*stringsize; i++)
		link[i] = 0;

	stateFileName = _stateFileName;

	allocated = true;
}

StateVariable::~StateVariable()
{
	statesize = 0;
	stringsize = 0;
	opCutOff = 0;
	nOp = 0;

	delete[] state;
	delete[] opString;
	delete[] vtx;
	delete[] link;
	offdiagOperatorLocs.clear();

	allocated = false;
}

void StateVariable::reset()
{
	if(allocated)
	{
		this->nOp = 0;

		for(long i = 0; i<stringsize; i++)
		{
			opString[i] = 0;
			vtx[i] = 0;
		}

		for(long i=0; i<4*stringsize; i++)
			link[i] = 0;
	}
}

void StateVariable::allocate()
{
	if(allocated)
	{
		delete[] state;
		delete[] opString;
		delete[] vtx;
		delete[] link;
	}

	state = new char[statesize];
	opString = new bigi[stringsize];
	vtx = new bigi[stringsize];
	link = new bigi[4*stringsize];

	reset();

	allocated = true;
}

int StateVariable::writeState(string fname)
{
	FILE* of = fopen(fname.c_str(),"wb");
	fwrite(&nOp,sizeof(nOp),1,of);
	fwrite(&opCutOff,sizeof(opCutOff),1,of);
	fwrite(&statesize,sizeof(statesize),1,of);
	fwrite(opString,sizeof(bigi),opCutOff,of);
	fwrite(state, sizeof(char),statesize,of);

	fwrite(&beta,sizeof(beta),1,of);
	fwrite(&lpContrIter,sizeof(lpContrIter),1,of);
	fclose(of);

	ibeta = 1.0/beta;

	return SUCCESS;
}

int StateVariable::reSize(long newSize)
{
	if(newSize>HARDCUTOFF)
	{
		fprintf(stdo,"Can't resize. Exceeds System Memory Requirements\n");
		return EXITCODE;
	}

	delete[] opString;
	delete[] vtx;
	delete[] link;
	stringsize = newSize;

	opString = new bigi[stringsize];
	vtx = new bigi[stringsize];
	link = new bigi[4*stringsize];

	for(long i = 0; i<stringsize; i++)
	{
		opString[i] = 0;
		vtx[i] = 0;
	}
	for(long i=0; i<4*stringsize; i++)
		link[i] = 0;

	return SUCCESS;
}

int StateVariable::readState(string fname)
{
	long localSize = 0;
	FILE* of = fopen(fname.c_str(),"rb");
	if(!of)
		return FILENOTFOUND;

	fread(&nOp,sizeof(nOp),1,of);
	fread(&localSize,sizeof(opCutOff),1,of);
	fread(&statesize,sizeof(statesize),1,of);

	if(localSize>stringsize)
	{
		fprintf(stdo,"Required Cutoff: %ld\n",localSize);
		fprintf(stdo,"Allocated size: %ld\n",stringsize);
		fprintf(stdo,"Exceeded Cutoff. Reallocating...\n");

		int status = reSize(localSize);
		if(status<0)
		{
			fprintf(stdo,"Exceeded System Memory Requirements\n");
			return EXITCODE;
		}
	}

	opCutOff = localSize;
	fread(opString,sizeof(bigi),opCutOff,of);
	fread(state, sizeof(char),statesize,of);
	fread(&beta,sizeof(beta),1,of);
	fread(&lpContrIter,sizeof(lpContrIter),1,of);
	fclose(of);

	ibeta = 1.0/beta;

	return SUCCESS;
}

void StateVariable::display(FILE* stdo)
{
	int nparticles = 0;
	for(int i=0;i<this->statesize;i++)
		nparticles += this->state[i];

	fprintf(stdo,"===================================================\n");
	fprintf(stdo,"State Variable Parameters\n");
	fprintf(stdo,"===================================================\n");
	fprintf(stdo,"File Name: %s\n",this->stateFileName.c_str());
	fprintf(stdo,"State Size: %d\n",this->statesize);
	fprintf(stdo,"Operator Cutoff: %ld\n",this->opCutOff);
	fprintf(stdo,"Maximum Allocated Size: %ld\n",this->stringsize);
	fprintf(stdo,"Number of Operators: %ld\n",this->nOp);
	fprintf(stdo,"Loop Control Iteration: %d\n",this->lpContrIter);
	fprintf(stdo,"Beta: %f\n",this->beta);
	fprintf(stdo,"No. of particles: %d\n",nparticles);
	fprintf(stdo,"\n\n");
	fflush(stdo);
}

int StateVariable::readState()
{
	return readState(this->stateFileName);
}

int StateVariable::writeState()
{
	return writeState(this->stateFileName);
}
