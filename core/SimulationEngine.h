/*
 * SimulationEngine.h
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

#ifndef SIMULATIONENGINE_H_
#define SIMULATIONENGINE_H_

namespace core
{
	class SimulationEngine
	{

	protected:
		long* first;
		long* last;
		char* state_t;

		long long maxLoopLength;
		long double invNumBonds;
		double multB;
		double imultB;

		gsl_rng *rgenref;

		Hamiltonian* hamiltonian;
		TransitionMap* tmap;
		StateVariable* stateVariable;

		FILE* stdo;

		///////////////////////////////////////////////////////
		friend class Hamiltonian;
		friend class StateVariable;
		///////////////////////////////////////////////////////

	public:

		SimulationEngine(Hamiltonian*, TransitionMap*, StateVariable*, gsl_rng*, FILE*);
		~SimulationEngine();


		void setRandomGenerator(gsl_rng* _lgen);
		void setLog(FILE* log);

		virtual int diagonalUpdate() = 0;
		virtual long loopUpdate(long long&) = 0;
		virtual int sweep(int&,long long&) = 0;
	};
}

#endif /* SIMULATIONENGINE_H_ */
