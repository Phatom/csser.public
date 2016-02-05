/*
 * Flags.h
 *
 *  Created on: Apr 25, 2014
 *      Author: ushnish

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

#ifndef FLAGS_H_
#define FLAGS_H_

const std::string runModeText[4] = {"Equilibriate","Pre-Warmup","Warmup","Data"};
const std::string boolean[2] = {"False","True"};

struct Flags
{
	bool oride_LoopContrIter;
	bool oride_temperature;

	float checkPointInterval;

	bool autoCalibrate;
	float autoFactor;
	int autoWindow;

	bool randomize;
	int seed;

	bool fixedGreenTime;
	float greendt;
	float sigma;
	bool outputTimeSeries;
	string timeSeriesFile;

	int loadFlags(string filename)
	{
		string temp;
		ifstream rif(filename);
		if(rif)
		{
			rif>>temp;
			if(temp.compare(boolean[0]) == 0)
				oride_LoopContrIter = false;
			else
				oride_LoopContrIter = true;

			rif>>temp;
			if(temp.compare(boolean[0]) == 0)
				oride_temperature = false;
			else
				oride_temperature = true;


			rif>>checkPointInterval;

			rif>>temp;
			if(temp.compare(boolean[1])==0)
				autoCalibrate = true;
			else
				autoCalibrate = false;
			rif>>autoFactor;
			rif>>autoWindow;

			rif>>temp;
			if(temp.compare(boolean[1])==0)
				randomize = true;
			else
				randomize = false;

			rif>>seed;

			rif>>temp;
			if(temp.compare(boolean[1])==0)
				fixedGreenTime = true;
			else
				fixedGreenTime = false;

			rif>>greendt;
			rif>>sigma;

			rif>>temp;
			if(temp.compare(boolean[1])==0)
				outputTimeSeries = true;
			else
				outputTimeSeries = false;

			rif>>timeSeriesFile;
		}
		else
			return EXITCODE;
		rif.close();

		return SUCCESS;
	}

	void display(FILE* stdo)
	{
		fprintf(stdo,"==========================================================================\n");
		fprintf(stdo,"Simulation Flags\n");
		fprintf(stdo,"==========================================================================\n");

		fprintf(stdo,"Override Loop Control: %s\n",boolean[(int) this->oride_LoopContrIter].c_str());
		fprintf(stdo,"Override Temperature: %s\n",boolean[(int) this->oride_temperature].c_str());
		fprintf(stdo,"Check Point Interval: %f\n",this->checkPointInterval);
		fprintf(stdo,"Auto Calibrate: %s\n",boolean[(int) this->autoCalibrate].c_str());
		fprintf(stdo,"Auto Calibrate factor: %f\n",this->autoFactor);
		fprintf(stdo,"Auto Window: %d\n",this->autoWindow);
		fprintf(stdo,"Randomize: %s\n",boolean[(int) this->randomize].c_str());
		fprintf(stdo,"Static Random Seed: %d\n",this->seed);
		fprintf(stdo,"Fix Green time step: %s\n",boolean[(int) this->fixedGreenTime].c_str());
		fprintf(stdo,"Green time step (delta_t): %f\n",this->greendt);
		fprintf(stdo,"Green sigma window: %f\n",this->sigma);
		fprintf(stdo,"Output Green time step stats: %s\n",boolean[(int) this->outputTimeSeries].c_str());
		fprintf(stdo,"Green time step stats. file: %s\n",timeSeriesFile.c_str());
		fprintf(stdo,"\n\n");
		fflush(stdo);
	}
};

#endif /* FLAGS_H_ */
