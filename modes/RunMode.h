/*
 * RunMode.h
 *
 *  Created on: Apr 23, 2014
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

#ifndef RUNMODE_H_
#define RUNMODE_H_

namespace runmode
{
	class RunMode
	{
	protected:
		core::Hamiltonian* hamiltonian;
		core::TransitionMap* tmap;
		core::StateVariable* runState;
		core::SimulationEngine* sengine;

		std::vector<measures::Measurable*> obsvCollection;

		string outputHeader;
	public:

		RunMode(core::Hamiltonian* h,core::TransitionMap* t,core::StateVariable* s,core::SimulationEngine* se)
		{
			hamiltonian = h;
			tmap = t;
			runState = s;
			sengine = se;
			outputHeader = "";
		}

		void setHamiltonian(core::Hamiltonian* h) {hamiltonian = h;}
		void setTransitionMap(core::TransitionMap* t) {tmap=t;}
		void setSimulationState(core::StateVariable* s) {runState = s;}
		void setSimulationEngine(core::SimulationEngine* eng) {sengine = eng;}
		void setOutputHeader(string _header) {outputHeader = _header;}
		std::vector<measures::Measurable*>& getObservableCollection() {return obsvCollection;}

		virtual int initialize() = 0;
		virtual void run() = 0;
		virtual void cleanup() = 0;
	};
}

#endif /* RUNMODE_H_ */
