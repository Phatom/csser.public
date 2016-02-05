/*
 * EquibRunMode_run.cpp
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

#include "mainincs.h"
#include "core.h"
#include "observables.h"
#include "RunMode.h"
#include "DataRunMode.h"

using namespace std;
using namespace __gnu_cxx;
using namespace runmode;
using namespace measures;

int DataRunMode::initialize()
{
	//Load equilibrium state and associated variables need for equilibration
	int ret = this->loadState();
	if(ret == EXITCODE)
	{
		fprintf(this->log,"Error in allocating memory\n");
		fflush(log);
		return EXITCODE;
	}
	else if(ret == FILENOTFOUND)
	{
		fprintf(log,"Data state found but corresponding core state file not found. Cannot do data gathering\n");
		fflush(log);
		return EXITCODE;
	}
	else if(ret == PARTIALFAIL)
	{
		fprintf(log,"Data state file not found. Starting new simulation.\n");
		fflush(log);
	}
	else
	{
		fprintf(log,"Data state found. Continuing new simulation.\n");
		for(int i=0;i<this->obsvCollection.size();i++)
			this->obsvCollection[i]->checkPointRead();
		fflush(log);
	}

	//Check flags over overrides
	if(flag->oride_LoopContrIter)
	{
		fprintf(this->log,"Overriding Loop Control from: %d to %d\n",runState->lpContrIter,this->params.startLoopControlIter);
		runState->lpContrIter = this->params.startLoopControlIter;
		fflush(this->log);
	}

	if(flag->oride_temperature)
	{
		fprintf(this->log,"Overriding temperature from: %f to %f\n",runState->beta,this->params.targetBeta);
		runState->beta = this->params.targetBeta;
		fflush(this->log);
	}

	initialized = true;
	return SUCCESS;
}

void DataRunMode::run()
{
	//Initialize here
	if(!initialized)
	{
		int ecode = this->initialize();
		if(ecode==EXITCODE)
			return;
	}

	fprintf(this->log,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
	this->runState->display(this->log);
	fprintf(this->log,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
	fflush(log);

	//Setup checkpointing
	int chkPInterval = int(flag->checkPointInterval*params.numSweeps + 0.5);
	chkPInterval = (chkPInterval == 0) ? 1 : chkPInterval;

	fprintf(log,"%s\n",this->outputHeader.c_str());
	//Start runs
	for(int cbin = this->currState.bin;cbin<this->params.numBins;cbin++)
	{
		for(int csweep=this->currState.sweep;csweep<this->params.numSweeps;csweep++)
		{
			fprintf(log,"%5d\t%5d\t",cbin,csweep);
			fflush(log);

			int der = 0;
			long long lup_er = 0;
			int error = this->sengine->sweep(der,lup_er);

			if(error==EXITCODE)
			{
				fprintf(this->log,"Error during sweep!\n");
				return;
			}

			///////////////////////////////////////////////////////////////////////////
			// After skipping do data gathering
			///////////////////////////////////////////////////////////////////////////

			eqbBase->measure();
			eqbBase->display();
			fflush(log);					

			for(int i=0;i<this->obsvCollection.size();i++)
				(this->obsvCollection)[i]->measure();

			///////////////////////////////////////////////////////////////////////////
			if(csweep==this->params.numSweeps-1 || csweep%chkPInterval==0 && csweep>0)
			{
				fprintf(log,"\tCKPT");
				fflush(log);

				currState.bin = cbin;
				currState.sweep = csweep+1;
				this->saveCheckPoint();

				for(int i=0;i<this->obsvCollection.size();i++)
					(this->obsvCollection)[i]->checkPointWrite();
			}

			fprintf(log,"\n");
			fflush(log);
			///////////////////////////////////////////////////////////////////////////

		}

		//Accumulate Data
		for(int i=0;i<this->obsvCollection.size();i++)
			(this->obsvCollection)[i]->writeViaIndex(cbin);

		//Now do calibration
		//1) Save temp file
		currState.bin = cbin + 1;
		currState.sweep = 0;
		this->saveState();

		//3) Reset basic observable
		eqbBase->clear();
	}

	this->cleanup();
	fprintf(log,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n");
	fflush(log);
}



