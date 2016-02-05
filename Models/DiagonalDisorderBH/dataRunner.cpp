/*
 * mainRunner.cpp
 *
 *  Created on: Apr 24, 2014
 *      Author: ushnish

 	Copyright (c) 2014 Ushnish Ray
	All rights reserved.

	This program and the accompanying materials are made available under explicit agreement
	between Ushnish Ray and the end user. You may not redistribute the code and
	accompanying materials to anyone.

	On the event that the software is used to generate data that is used implicitly or explicitly
	for research purposes, proper acknowledgment must be provided  in the citations section of
	publications.

	This software cannot be used for commercial purposes in any way whatsoever.
 */

#include "sse.h"
#include "StateVariableBH.h"
#include "BoseHubbardDDisorder.h"
#include "EqbSimEngine.h"
#include "BHspdmSimEngine.h"
#include "BHGreenSimEngine.h"

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

	string dataParamFile;
	string staticObsFile;
	string derivedObsFile;
	string localObsFile;
	string greenFile;
	string sfluidFile;
	string sfluidLocalFile;

	void display(ofstream &outl)
	{
		outl<<"------------------------------------------------------------\n";
		outl<<"Data Session Files\n";
		outl<<"------------------------------------------------------------\n";
		outl<<"Output File: "<<outputFile<<endl;
		outl<<"Flag File: "<<flagFile<<endl;
		outl<<"State File: "<<startStateFile<<endl;
		outl<<"Hamiltonian File: "<<hamParamFile<<endl;
		outl<<"Base Directory: "<<baseDir<<endl;
		outl<<"Session Parameter File: "<<dataParamFile<<endl;
		outl<<endl;
		outl<<"Observable File Names:\n";
		outl<<"Static Observable: "<<staticObsFile<<endl;
		outl<<"Derived Observable: "<<derivedObsFile<<endl;
		outl<<"Local Observable: "<<localObsFile<<endl;
		outl<<"Green Observable: "<<greenFile<<endl;
		outl<<"Superfluid Observable: "<<sfluidFile<<endl;
		outl<<"Superfluid Local Observable: "<<sfluidLocalFile<<endl;
		outl<<"------------------------------------------------------------\n";
	}
};

int dataRunner(int rank, string path, char* argv)
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
	SimulationEngine* simengine;

	////////////////////////////////////////////////
	fnames.mainParamFile = path + string(argv);
	ifstream rif(fnames.mainParamFile);
	if(rif)
	{
		rif>>fnames.baseDir; fnames.baseDir = path + fnames.baseDir;
		rif>>fnames.outputFile; fnames.outputFile = path + fnames.outputFile;
		rif>>fnames.flagFile; fnames.flagFile = path + fnames.flagFile;
		rif>>fnames.startStateFile; fnames.startStateFile = path + fnames.startStateFile;
		rif>>fnames.hamParamFile; fnames.hamParamFile = path + fnames.hamParamFile;
		rif>>fnames.stateOutputFile; fnames.stateOutputFile = path + fnames.stateOutputFile;

		rif>>fnames.dataParamFile; fnames.dataParamFile = path + fnames.dataParamFile;
		rif>>fnames.staticObsFile; fnames.staticObsFile = path + fnames.staticObsFile;
		rif>>fnames.derivedObsFile; fnames.derivedObsFile = path + fnames.derivedObsFile;
		rif>>fnames.localObsFile; fnames.localObsFile = path + fnames.localObsFile;
		rif>>fnames.greenFile; fnames.greenFile = path + fnames.greenFile;
		rif>>fnames.sfluidFile; fnames.sfluidFile = path + fnames.sfluidFile;
		rif>>fnames.sfluidLocalFile; fnames.sfluidLocalFile = path + fnames.sfluidLocalFile;
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
	simflags.timeSeriesFile = path + simflags.timeSeriesFile;
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
	try{
		svar = new StateVariableBH(stdo,fnames.startStateFile,fnames.stateOutputFile);
	}
	catch(int exception)
	{
		if(exception == ConstructionException)
		{
			fprintf(stdo,"String Size is larger than HARDCUTOFF. Re-compile after increasing size.\n");
			fflush(stdo);
			return EXITCODE;
		}
	}
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

	//6) Setup Simulation Engine based on mode desired
	DataParameters dataparam;
	error = dataparam.loadDataParams(fnames.dataParamFile);
	if(error==EXITCODE)
	{
		fprintf(stdo,"Error load run parameters file.\n");
		fclose(stdo);
		delete svar;
		delete bhamiltonian;
		delete simengine;
		return EXITCODE;
	}
	dataparam.dataStateFile = path + dataparam.dataStateFile;
	dataparam.display(stdo);

	//7) Setup observables
	map2d spdm;
	GijI lpvargreen;
	if(dataparam.mode == NONE)
	{
		simengine = new EqbSimEngine(bhamiltonian,tmap,svar,rgenref,stdo);
	}
	else if(dataparam.mode == SPDM || dataparam.mode == SPDMANDSFLUID)
	{
		simengine = new BHspdmSimEngine(bhamiltonian,tmap,svar,rgenref,stdo,spdm);
	}
	else if(dataparam.mode == GREEN || dataparam.mode == ALL)
	{
		simengine = new BHGreenSimEngine(bhamiltonian,tmap,svar,rgenref,stdo,lpvargreen);
	}

	//8) Setup Runners
	DataRunMode* datarun = new DataRunMode(bhamiltonian,tmap,svar,simengine,dataparam,&simflags,stdo);
	datarun->setOutputHeader("Bin\t\tStep\t\tDiag. Time\t\tLoop Time\t\tOperator Count\t\tDensity\t\tEnergy\n");

	int ecode = datarun->initialize(); //This is really important for many state dependent observables like green's function
	if(ecode==EXITCODE)
		return EXITCODE;

	//9) Setup Observables
	std::vector<measures::Measurable*> &obsvCollection = datarun->getObservableCollection();
	BasicGlobalObservable globalData(fnames.staticObsFile,svar,bhamiltonian,stdo);
	DerivedGlobalObservable globalDerived(fnames.derivedObsFile,svar,bhamiltonian,stdo);
	BasicLocalObservable localData(fnames.localObsFile,svar,bhamiltonian,stdo);


	//Super-fluid observables are special
	//Setup state list
	SfluidWindingObservable sfluidData(fnames.sfluidFile,svar,bhamiltonian,stdo);
	SfluidRotationalObservable sfluidRData(fnames.sfluidFile,svar,bhamiltonian,stdo);
	SPDMObservable spdmData(fnames.greenFile,svar,bhamiltonian,stdo,spdm);
	GreenObservable greenData(fnames.greenFile,svar,bhamiltonian,stdo,simflags,lpvargreen);

	obsvCollection.push_back(&globalData);
	obsvCollection.push_back(&globalDerived);
	obsvCollection.push_back(&localData);

	if(dataparam.mode == SPDM)
	{
		obsvCollection.push_back(&spdmData);
	}
	else if(dataparam.mode == SPDMANDSFLUID)
	{
		if(bhamiltonian->bc == HW)
			obsvCollection.push_back(&sfluidRData);
		else
			obsvCollection.push_back(&sfluidData);

		obsvCollection.push_back(&spdmData);
	}
	else if(dataparam.mode == GREEN)
	{
		greenData.createTvalues();
		obsvCollection.push_back(&greenData);
	}
	else if(dataparam.mode == ALL)
	{
		if(bhamiltonian->bc == HW)
			obsvCollection.push_back(&sfluidRData);
		else
			obsvCollection.push_back(&sfluidData);

		greenData.createTvalues();
		obsvCollection.push_back(&greenData);
	}

	fprintf(stdo,"Observables setup successfully\n");
	fflush(stdo);
	/////////////////////////////////////////////////////
	// Run Modes
	/////////////////////////////////////////////////////
	datarun->run();

	/////////////////////////////////////////////////////
	// Cleanup
	/////////////////////////////////////////////////////
	fprintf(stdo,"Complete.\n\n");
	fclose(stdo);

	delete svar;
	delete bhamiltonian;
	delete tmap;
	delete simengine;
	delete datarun;
}

