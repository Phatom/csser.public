/*
 * mainRunner.cpp
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

#include "sse.h"
#include "StateVariableBH.h"
#include "BoseHubbardDDisorder.h"
#include "EqbSimEngine.h"

using namespace std;
using namespace core;
using namespace measures;
using namespace runmode;

struct FileNames
{
	string mainParamFile;

	string baseDir;
	string outputFile;
	string flagFile;
	string startStateFile;
	string hamParamFile;
	string stateOutputFile;

	string runParamFile;

	void display(ofstream &llog)
	{
		llog<<"------------------------------------------------------------\n";
		llog<<"Warm Session Files\n";
		llog<<"------------------------------------------------------------\n";
		llog<<"Output File: "<<outputFile<<endl;
		llog<<"Flag File: "<<flagFile<<endl;
		llog<<"State File: "<<startStateFile<<endl;
		llog<<"Hamiltonian File: "<<hamParamFile<<endl;
		llog<<"Base Directory: "<<baseDir<<endl;
		llog<<"Session Parameter File: "<<runParamFile<<endl;
		llog<<"------------------------------------------------------------\n";
	}
};

int warmRunner(int rank, string path, char* argv)
{
	FileNames fnames;

	////////////////////////////////////////////////
	//Setup Random generator
	gsl_rng_env_setup();
	gsl_rng* rgenref = gsl_rng_alloc(gsl_rng_mt19937);
	////////////////////////////////////////////////

	Flags simflags;
	StateVariableBH* svar;
	BoseHubbardDDisorder* bhamiltonian;
	SSTransitionMap* tmap;
	EqbSimEngine*	simengine;

	////////////////////////////////////////////////
	fnames.mainParamFile = string(argv);
	ifstream rif(argv);
	if(rif)
	{
		rif>>fnames.baseDir; fnames.baseDir = path + fnames.baseDir;
		rif>>fnames.outputFile; fnames.outputFile = path + fnames.outputFile;
		rif>>fnames.flagFile; fnames.flagFile = path + fnames.flagFile;
		rif>>fnames.startStateFile; fnames.startStateFile = path + fnames.startStateFile;
		rif>>fnames.hamParamFile; fnames.hamParamFile = path + fnames.hamParamFile;
		rif>>fnames.stateOutputFile; fnames.stateOutputFile = path + fnames.stateOutputFile;

		rif>>fnames.runParamFile; fnames.runParamFile = path + fnames.runParamFile;
	}
	else
	{
		cout<<fnames.mainParamFile<<" not found!\n";
		return EXITCODE;
	}
	rif.close();

	ofstream llog(fnames.outputFile,ios::app);
	fnames.display(llog);
	llog.close();

	////////////////////////////////////////////////

	/////////////////////////////////////////////////////
	// Setup
	/////////////////////////////////////////////////////

	//1) Process Flags and files
	FILE* stdo = fopen(fnames.outputFile.c_str(),"a");
	int error = simflags.loadFlags(fnames.flagFile);
	if(error==EXITCODE)
	{
		fprintf(stdo,"Error loading flag file.\n");
		fclose(stdo);
		return EXITCODE;
	}
	simflags.display(stdo);

	//2) Setup Random Number Generators
	if(simflags.randomize)
		simflags.seed = (int) time(NULL) + rank;
	else
		simflags.seed += rank;
	fprintf(stdo,"Seed: %d\n",simflags.seed);

	gsl_rng_set(rgenref,simflags.seed);
	fprintf(stdo,"Test: %d\n",gsl_rng_uniform_int(rgenref,100));
	fprintf(stdo,"Test: %d\n",gsl_rng_uniform_int(rgenref,100));
	fflush(stdo);

	//3) Setup State Variable
	svar = new StateVariableBH(stdo,fnames.startStateFile,fnames.stateOutputFile);
	svar->rgenref = rgenref;
	svar->initializeState();
	svar->display(stdo);

	//4) Setup Hamiltonian
	bhamiltonian = new BoseHubbardDDisorder(rank,path,fnames.hamParamFile,stdo);
	error = bhamiltonian->loadHamiltonianParams();
	if(error==EXITCODE)
	{
		fprintf(stdo,"Error hamiltonian parameters file.\n");
		fclose(stdo);
		delete svar;
		delete bhamiltonian;
		return EXITCODE;
	}
	bhamiltonian->initialize();
	bhamiltonian->display();

	//5) Setup Transition File
	tmap = new SSTransitionMap(stdo,bhamiltonian);
	tmap->initialize(fnames.baseDir);

	//6) Setup Simulation Engine
	simengine = new EqbSimEngine(bhamiltonian,tmap,svar,rgenref,stdo);

	//7) Setup Runners
	WarmParameters warmparam;
	error = warmparam.loadWarmParams(fnames.runParamFile);
	if(error==EXITCODE)
	{
		fprintf(stdo,"Error load run parameters file.\n");
		fclose(stdo);
		delete svar;
		delete bhamiltonian;
		delete simengine;
		return EXITCODE;
	}
	warmparam.warmStateFile = path + warmparam.warmStateFile;
	warmparam.display(stdo);

	WarmRunMode* warmrun = new WarmRunMode(bhamiltonian,tmap,svar,simengine,warmparam,&simflags,stdo);
	warmrun->setOutputHeader("Step\t\tVertex Count\t\tBounce Count\t\tLoops Used\t\tOperator Count\t\tDensity\t\tEnergy\n");
	/////////////////////////////////////////////////////
	// Run Modes
	/////////////////////////////////////////////////////
	warmrun->run();

	/////////////////////////////////////////////////////
	// Cleanup
	/////////////////////////////////////////////////////
	fprintf(stdo,"Complete.\n\n");
	fclose(stdo);

	delete svar;
	delete bhamiltonian;
	delete tmap;
	delete simengine;
	delete warmrun;
}

