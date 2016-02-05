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

#ifndef EQUIBRUNMODE_H_
#define EQUIBRUNMODE_H_

namespace runmode
{
	struct EqbParameters
	{
		float startBeta;
		float incBeta;
		float targetBeta;

		int numSweeps;
		int skipSweeps;
		int startLoopControlIter;

		string eqbStateFile;

		int loadEqbParams(std::string fn)
		{
			ifstream rif(fn.c_str());
			if(rif)
			{
				rif>>startBeta;
				rif>>incBeta;
				rif>>targetBeta;
				rif>>numSweeps;
				rif>>skipSweeps;
				rif>>startLoopControlIter;
				rif>>eqbStateFile;
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
			fprintf(stdo,"Start Beta: %f\n",startBeta);
			fprintf(stdo,"Inc Beta: %f\n",incBeta);
			fprintf(stdo,"Target Beta: %f\n",targetBeta);
			fprintf(stdo,"Number of Sweeps: %d\n",numSweeps);
			fprintf(stdo,"Number of Sweeps to skip: %d\n",skipSweeps);
			fprintf(stdo,"Loop Control Iteration: %d\n",startLoopControlIter);
			fprintf(stdo,"\n\n");
			fflush(stdo);
		}
	};

	struct EqbState
	{
		float beta;
		int sweep;
	};

	class EquibRunMode : public RunMode
	{
	protected:
		FILE* log;

		EqbParameters eparams;
		EqbState currState;
		Flags* flag;

		measures::BasicGlobalObservable* eqbBase;
	public:

		EquibRunMode(core::Hamiltonian*,core::TransitionMap*,core::StateVariable*,core::SimulationEngine*, EqbParameters& , Flags*, FILE*);
		~EquibRunMode();

		void setParameters(EqbParameters& ep) {eparams = ep;}
		void setFlags(Flags* f) {flag = f;}

		int loadEqbState();
		void saveEqbState();
		void saveCheckPoint();
		void saveTemperatureState(int);

		int initialize();
		void run();
		void cleanup() {};
	};
}

#endif /* EQUIBRUNMODE_H_ */
