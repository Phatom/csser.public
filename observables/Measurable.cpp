/*
 * Measurable.cpp
 *
 *  Created on: Apr 22, 2014
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

Measurable::Measurable(string fn, core::StateVariable* lstate, core::Hamiltonian* lham, FILE* llog)
{
	fileName = fn;
	chkFileName = fileName + ".chk";
	state = lstate;
	hamiltonian = lham;
	log = llog;
}

void Measurable::setLog(FILE* llog)
{
	log = llog;
}

void Measurable::setState(core::StateVariable *s)
{
	state = s;
}

void Measurable::setFileName(string fname)
{
	fileName = fname;
}
