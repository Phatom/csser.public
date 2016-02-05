/*
 * StateVariableEQB.h
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

#ifndef STATEVARIABLEEQB_H_
#define STATEVARIABLEEQB_H_

#include "StateVariable.h"

class StateVariableBH: public core::StateVariable {
public:

	gsl_rng* rgenref;
	int LAM;

	StateVariableBH(FILE* _stdo,string startStateFile,string saveLocation);

	int initializeState();
};

#endif /* STATEVARIABLEEQB_H_ */
