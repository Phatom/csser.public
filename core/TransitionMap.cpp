/*
 * TranstionMap.cpp
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

using namespace std;
using namespace __gnu_cxx;
using namespace core;

void TransitionMap::writeBondProbs(string base)
{
	  FILE* tp=fopen((base + "TPROBS").c_str(),"wb");
	 fwrite(bondToProb,sizeof(double),tsize,tp);
	 fclose(tp);
}

int TransitionMap::readBondProbs(string base)
{
	 FILE* tp=fopen((base + "TPROBS").c_str(),"rb");
	 if(!tp)
		 return FILENOTFOUND;

	 fread(bondToProb,sizeof(double),tsize,tp);
	 fclose(tp);

	 //displayBondProbs();
	 return SUCCESS;
}

void TransitionMap::initialize(string base)
{
	this->initializeMaps();
	this->calcExitProbs(base);
}
