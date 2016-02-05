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
#include "WarmRunMode.h"

using namespace std;
using namespace __gnu_cxx;
using namespace runmode;
using namespace measures;

int WarmRunMode::initialize()
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
		fprintf(log,"Warm state found but corresponding core state file not found. Cannot do warmup\n");
		fflush(log);
		return EXITCODE;
	}
	else if(ret == PARTIALFAIL)
	{
		fprintf(log,"Warm state file not found. Starting new simulation.\n");
		fflush(log);
	}
	else
	{
		fprintf(log,"Warm states found. Continuing new simulation.\n");

		for(int i=0;i<this->obsvCollection.size();i++)
			this->obsvCollection[i]->checkPointRead();

		fflush(log);
	}
	//Check flags over overrides
	if(flag->oride_LoopContrIter)
		runState->lpContrIter = this->params.startLoopControlIter;

	if(flag->oride_temperature)
		runState->beta = this->params.targetBeta;

	return SUCCESS;
}

void WarmRunMode::run()
{
	//Initialize here
	int ecode = this->initialize();
	if(ecode==EXITCODE)
		return;

	//Setup checkpointing
	int chkPInterval = int(flag->checkPointInterval*params.numSweeps + 0.5);

	//Start runs
	fprintf(log,"\n");
	fprintf(this->log,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
	fprintf(this->log,"Current Temperature Sweep\n");
	fprintf(this->log,"Beta: %f\n",this->runState->beta);
	this->runState->display(this->log);
	fprintf(this->log,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
	fflush(log);

	long long vtxVisitedCount = 0, vtxVisitedCount2 = 0;
	int npts = 0;
	long autovtxCount = 0;
	int autonpts = 0;

	fprintf(log,"%s\n",this->outputHeader.c_str());
	for(int csweep=this->currState.sweep;csweep<this->params.numSweeps;csweep++)
	{
		fprintf(log,"%5d\t",csweep);

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
		vtxVisitedCount += lup_er;
		vtxVisitedCount2 += lup_er*lup_er;
		npts++;

		if(flag->autoCalibrate)
		{
			autovtxCount += lup_er;
			autonpts++;
		}

		eqbBase->measure();
		eqbBase->display();


		for(int i=0;i<this->obsvCollection.size();i++)
			(this->obsvCollection)[i]->measure();

		///////////////////////////////////////////////////////////////////////////
		if(csweep==this->params.numSweeps-1 || csweep%chkPInterval==0 && csweep>0)
		{
			fprintf(log,"\tCKPT");
			currState.sweep = csweep+1;
			this->saveCheckPoint();

			for(int i=0;i<this->obsvCollection.size();i++)
				(this->obsvCollection)[i]->checkPointWrite();
		}

		fprintf(log,"\n");
		fflush(log);
		///////////////////////////////////////////////////////////////////////////
		if(flag->autoCalibrate && autonpts==flag->autoWindow)
		{
			long cmpWith = flag->autoFactor*eqbBase->staticO.operatorCount;
			double avgVtxCount = double(autovtxCount)/autonpts;

			if(cmpWith - avgVtxCount >  0.1*cmpWith) //to within 10% is OK otherwise adjust
			{
				if(runState->lpContrIter > 10)
					runState->lpContrIter *= 1.15; //increase by some proportion
				else
					runState->lpContrIter++;
			}
			else if(runState->lpContrIter>1 && avgVtxCount - cmpWith > 0.1*cmpWith)
			{
				runState->lpContrIter *= 0.9; //Reduce to 90%
				if(runState->lpContrIter == 0)
					runState->lpContrIter = 1;
			}

			autonpts = 0;
			autovtxCount = 0;
		}
	}

	//Accumulate Data
	for(int i=0;i<this->obsvCollection.size();i++)
		(this->obsvCollection)[i]->write();

	//1) Save temp file
	this->saveState();

	//Calculate averages
	double avgVtxCount = double(vtxVisitedCount)/npts, avgVtxCount2 = sqrt((long double)(vtxVisitedCount2)/npts - avgVtxCount*avgVtxCount);
	double avg_nop = 0.0, avg_rho = 0.0, avg_energy = 0.0;
	double avg_nop2 = 0.0, avg_rho2 = 0.0, avg_energy2 = 0.0;
	int sample = 0;
	for(vector<staticObservable>::iterator it=eqbBase->collection.begin(); it != eqbBase->collection.end(); ++it)
	{
		avg_nop += it->operatorCount; avg_nop2 += (long) it->operatorCount*it->operatorCount;
		avg_rho += it->density; avg_rho2 += it->density*it->density;
		avg_energy += it->energy; avg_energy2 += it->energy*it->energy;
		sample++;
	}
	avg_nop /= sample; avg_nop2 /= sample; avg_nop2 = sqrt(avg_nop2-avg_nop*avg_nop);
	avg_rho /= sample; avg_rho2 /= sample; avg_rho2 = sqrt(avg_rho2-avg_rho*avg_rho);
	avg_energy /= sample; avg_energy2 /= sample; avg_energy2 = sqrt(avg_energy2-avg_energy*avg_energy);

	fprintf(this->log,"===========================================================\n");
	fprintf(this->log,"Warmup Done\n");
	fprintf(this->log,"Number of Operators [Variance^0.5]: %16.10le [ %16.10le ]\n",avg_nop,avg_nop2);
	fprintf(this->log,"Average Vertex Count [Variance^0.5]: %16.10le [%16.10le ]\n",avgVtxCount,avgVtxCount2);
	fprintf(this->log,"Density [Variance^0.5]: %16.10le [ %16.10le ]\n",avg_rho,avg_rho2);
	fprintf(this->log,"Energy [Variance^0.5]: %16.10le [ %16.10le ]\n",avg_energy,avg_energy2);
	fprintf(this->log,"Current Loop Control Number: %d\n",this->runState->lpContrIter);

	fprintf(this->log,"===========================================================\n");
	fflush(log);

	this->cleanup();
	fprintf(log,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n");
	fflush(log);
}



