/*
 * basicSimEngine.h
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

#ifndef BASICSIMENGINE_H_
#define BASICSIMENGINE_H_

namespace core
{
	class BasicSimEngine:public SimulationEngine
	{

	public:

		BasicSimEngine(Hamiltonian* h, TransitionMap* tmap, StateVariable* s, gsl_rng* lgsl,
				FILE* llog):SimulationEngine(h,tmap,s,lgsl,llog) { };

		~BasicSimEngine() { }

		int diagonalUpdate();
		long loopUpdate(long long&);
	};
}

#endif /* BASICSIMENGINE_H_ */
