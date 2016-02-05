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

#ifndef DATARUNMODE_H_
#define DATARUNMODE_H_

namespace runmode
{
	enum DataGatherType {NONE,SPDM,SFLUID,GREEN,SPDMANDSFLUID,ALL};
	const std::string dataGatherModeText[6] = {"none","spdm","sfluid","green","spdm&sfluid","all"};

	struct DataParameters
	{
		float targetBeta;
		int numBins;
		int numSweeps;
		int startLoopControlIter;
		string dataStateFile;
		DataGatherType mode;

		int loadDataParams(std::string fn)
		{
			string smode;

			ifstream rif(fn.c_str());
			if(rif)
			{
				rif>>targetBeta;
				rif>>numBins;
				rif>>numSweeps;
				rif>>startLoopControlIter;
				rif>>dataStateFile;
				rif>>smode;
				rif.close();
			}
			else
				return EXITCODE;

			for(int i=0;i<6;i++)
				if(smode.compare(dataGatherModeText[i]) == 0)
					mode = static_cast<DataGatherType>(i);

			return SUCCESS;
		}

		void display(FILE* stdo)
		{
			fprintf(stdo,"===================================================\n");
			fprintf(stdo,"Data Parameters\n");
			fprintf(stdo,"===================================================\n");
			fprintf(stdo,"Target Beta: %f\n",targetBeta);
			fprintf(stdo,"Number of Bins: %d\n",numBins);
			fprintf(stdo,"Number of Sweeps: %d\n",numSweeps);
			fprintf(stdo,"Loop Control Iteration: %d\n",startLoopControlIter);
			fprintf(stdo,"Data state file: %s\n",dataStateFile.c_str());
			fprintf(stdo,"Data gathering mode: %s\n",dataGatherModeText[(int) mode].c_str());
			fprintf(stdo,"\n\n");
			fflush(stdo);
		}
	};

	struct DataState
	{
		int bin;
		int sweep;
	};

	class DataRunMode : public RunMode
	{
	protected:
		bool initialized;

		FILE* log;

		DataParameters params;
		DataState currState;
		Flags* flag;

		measures::BasicStaticObservable* eqbBase;
	public:

		DataRunMode(core::Hamiltonian*,core::TransitionMap*,core::StateVariable*,core::SimulationEngine*, DataParameters& , Flags*, FILE*);
		~DataRunMode();

		void setParameters(DataParameters& ep) {params = ep;}
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
