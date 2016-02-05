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
#include "EquibRunMode.h"

using namespace std;
using namespace __gnu_cxx;
using namespace runmode;
using namespace measures;

int EquibRunMode::initialize()
{
	//Load equilibrium state and associated variables need for equilibration
	int ret = this->loadEqbState();
	if(ret == EXITCODE)
	{
		fprintf(this->log,"Error in allocating memory\n");
		fflush(log);
		return EXITCODE;
	}
	else if(ret == FILENOTFOUND)
	{
		runState->lpContrIter = this->eparams.startLoopControlIter;
		fprintf(log,"Old state file not found. Starting new simulation\n");
		fflush(log);
	}
	else
	{
		fprintf(log,"Old state file found. Continuing simulation\n");

		for(int i=0;i<this->obsvCollection.size();i++)
			this->obsvCollection[i]->checkPointRead();

		fflush(log);
	}

	//Check flags over overrides
	if(flag->oride_LoopContrIter)
	{
		fprintf(this->log,"Overriding Loop Control from: %d to %d\n",runState->lpContrIter,eparams.startLoopControlIter);
		runState->lpContrIter = this->eparams.startLoopControlIter;
		fflush(this->log);
	}

	return SUCCESS;
}

void EquibRunMode::run()
{
	//Initialize here
	int ecode = this->initialize();
	if(ecode==EXITCODE)
		return;

	//Setup checkpointing
	int chkPInterval = int(flag->checkPointInterval*eparams.numSweeps + 0.5);


	//Check if last state was completed state
	if(this->currState.sweep>=this->eparams.numSweeps)
	{
		currState.sweep = 0;
		if(flag->autoCalibrate)
		{
			double bratio = 1.0+eparams.incBeta/runState->beta;
			
			long expectedNop = runState->nOp*bratio;
			if(expectedNop > runState->opCutOff)
				runState->opCutOff = 1.25*runState->opCutOff; //Increase it to accommodate next higher temperature
			else if(1.25*expectedNop < runState->opCutOff)
				runState->opCutOff = 1.25*expectedNop;

			/* BE VERY CAREFUL HERE
			 * Generally if we truncate the operator string (and only if we truncate)
			 * and then do not reset the string to 0. Then it is possible that some
			 * off-diagonal terms will be truncated that are essential to maintain pbc in time.
			 * SO if we truncate then we MUST reset. On the other hand if we do not truncate
			 * then there is no need to reset. Generally for equilibration purposes IT IS A GOOD IDEA
			 * to reset the string. Keep the state however.
			 */

#ifdef BOUNDCHECK
			if(runState->opCutOff>runState->stringsize)
				runState->reSize(runState->opCutOff);
#endif

			runState->lpContrIter *= 0.5/(bratio*bratio);
			if(runState->lpContrIter == 0)
				runState->lpContrIter = 1;
		}
	
		this->currState.beta += this->eparams.incBeta;	
		runState->reset();
		eqbBase->clear();		
	}

	//Start runs
	for(float cbeta = this->currState.beta;cbeta<=this->eparams.targetBeta;cbeta+=this->eparams.incBeta)
	{
		this->runState->beta = cbeta;
		this->runState->ibeta = 1.0/cbeta;

		currState.beta = cbeta;

		fprintf(log,"\n");
		fprintf(this->log,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
		fprintf(this->log,"Current Temperature Sweep\n");
		fprintf(this->log,"Beta: %f\n",cbeta);
		this->runState->display(this->log);
		fprintf(this->log,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
		fflush(log);

		//Set current temperature etc. for simulation and save new run
		long long vtxVisitedCount = 0, vtxVisitedCount2 = 0;
		int npts = 0;
		long autovtxCount = 0;
		int autonpts = 0;

		fprintf(log,"%s\n",this->outputHeader.c_str());
		for(int csweep=this->currState.sweep;csweep<this->eparams.numSweeps;csweep++)
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

			if(csweep<eparams.skipSweeps)
			{
				fprintf(log,"\n");
				continue;
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

			/*
			for(int i=0;i<this->obsvCollection.size();i++)
				(this->obsvCollection)[i].measure();
			*/

			///////////////////////////////////////////////////////////////////////////
			if(csweep==this->eparams.numSweeps-1 || csweep%chkPInterval==0 && csweep>0)
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

		int tmpoffset = floor(cbeta/eparams.incBeta);

		//Accumulate Data
		for(int i=0;i<this->obsvCollection.size();i++)
			(this->obsvCollection)[i]->writeViaIndex(tmpoffset);

		//Now do calibration
		//1) Save temp file
		this->saveEqbState();
		this->saveTemperatureState(tmpoffset);

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
		fprintf(this->log,"Current Temperature Done\n");
		fprintf(this->log,"Number of Operators [Variance^0.5]: %16.10le [ %16.10le ]\n",avg_nop,avg_nop2);
		fprintf(this->log,"Average Vertex Count [Variance^0.5]: %16.10le [%16.10le ]\n",avgVtxCount,avgVtxCount2);
		fprintf(this->log,"Density [Variance^0.5]: %16.10le [ %16.10le ]\n",avg_rho,avg_rho2);
		fprintf(this->log,"Energy [Variance^0.5]: %16.10le [ %16.10le ]\n",avg_energy,avg_energy2);
		fprintf(this->log,"Current Loop Control Number: %d\n",this->runState->lpContrIter);

		fprintf(this->log,"===========================================================\n");
		fflush(log);
		
		//Check if done.
		if(fabs(cbeta-eparams.targetBeta)<1.0e-6)
			break;

		//2) Calibrate
		if(flag->autoCalibrate)
		{
			double bratio = 1.0+eparams.incBeta/cbeta;
			long expectedNop = avg_nop*bratio;
			if(expectedNop > runState->opCutOff)
				runState->opCutOff = 1.25*runState->opCutOff; //Increase it to accommodate next higher temperature
			else if(1.25*expectedNop < runState->opCutOff)
				runState->opCutOff = 1.25*expectedNop;

			/* BE VERY CAREFUL HERE
			 * Generally if we truncate the operator string (and only if we truncate)
			 * and then do not reset the string to 0. Then it is possible that some
			 * off-diagonal terms will be truncated that are essential to maintain pbc in time.
			 * SO if we truncate then we MUST reset. On the other hand if we do not truncate
			 * then there is no need to reset. Generally for equilibration purposes IT IS A GOOD IDEA
			 * to reset the string. Keep the state however.
			 */

#ifdef BOUNDCHECK
			if(runState->opCutOff>runState->stringsize)
				runState->reSize(runState->opCutOff);
#endif

			runState->lpContrIter *= 0.5/(bratio*bratio);
			if(runState->lpContrIter == 0)
				runState->lpContrIter = 1;
		}

		//3) Finally Reset
		currState.sweep = 0;
		runState->reset();
		eqbBase->clear();
	}

	this->cleanup();
	fprintf(log,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n");
	fflush(log);
}



