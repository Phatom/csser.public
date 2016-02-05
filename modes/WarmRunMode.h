/*
 * EquibRunMode.h
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

#ifndef WARMRUNMODE_H_
#define WARMRUNMODE_H_

namespace runmode
{
	struct WarmParameters
	{
		float targetBeta;
		int numSweeps;
		int startLoopControlIter;
		string warmStateFile;

		int loadWarmParams(std::string fn)
		{
			ifstream rif(fn.c_str());
			if(rif)
			{
				rif>>targetBeta;
				rif>>numSweeps;
				rif>>startLoopControlIter;
				rif>>warmStateFile;
				rif.close();
			}
			else
				return EXITCODE;

			return SUCCESS;
		}

		void display(FILE* stdo)
		{
			fprintf(stdo,"===================================================\n");
			fprintf(stdo,"Equlibration Parameters\n");
			fprintf(stdo,"===================================================\n");
			fprintf(stdo,"Target Beta: %f\n",targetBeta);
			fprintf(stdo,"Number of Sweeps: %d\n",numSweeps);
			fprintf(stdo,"Loop Control Iteration: %d\n",startLoopControlIter);
			fprintf(stdo,"\n\n");
			fflush(stdo);
		}
	};

	struct WarmState
	{
		int sweep;
	};

	class WarmRunMode : public RunMode
	{
	protected:
		FILE* log;

		WarmParameters params;
		WarmState currState;
		Flags* flag;

		measures::BasicGlobalObservable* eqbBase;
	public:

		WarmRunMode(core::Hamiltonian*,core::TransitionMap*,core::StateVariable*,core::SimulationEngine*, WarmParameters& , Flags*, FILE*);
		~WarmRunMode();

		void setParameters(WarmParameters& ep) {params = ep;}
		void setFlags(Flags* f) {flag = f;}

		int loadState();
		void saveState();
		void saveCheckPoint();

		int initialize();
		void run();
		void cleanup() {};
	};
}

#endif /* EQUIBRUNMODE_H_ */
