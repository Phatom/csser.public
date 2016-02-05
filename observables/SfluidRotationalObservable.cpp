/*
 * SfluidRotationalObservable.cpp
 *
 *  Created on: Aug 7, 2014
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

namespace measures {

SfluidRotationalObservable::SfluidRotationalObservable(string s, core::StateVariable* lstate,
core::Hamiltonian* lham, FILE* llog) : Measurable(s,lstate,lham,llog)
{
	size = lham->numSites;

	//Stores int -> (z,y,x)
	int C = lham->nspd/2;
	siteMap = new int[3*size];
	int L2 = lham->nspd*lham->nspd;
	for(int i=0;i<size;i++)
	{
		int z = i/L2;
		int y = i%L2/lham->nspd;
		int x = i%L2%lham->nspd;

		siteMap[3*i] = z-C;
		siteMap[3*i+1] = y-C;
		siteMap[3*i+2] = x-C;
	}

	rho_s = new double[6*size];
	Atotal_r = new double[3*size];
	Atotal = new double[3];

	for(int i=0;i<3*size;i++)
		Atotal_r[i] = 0.0;

	for(int i=0;i<6*size;i++)
		rho_s[i] = 0.0;

	Zcount = 0;
	Atotal[0] = Atotal[1] = Atotal[2] = 0.0;
}

SfluidRotationalObservable::~SfluidRotationalObservable()
{
	size = 0;
	delete[] rho_s;
	delete[] Atotal;
	delete[] Atotal_r;
	delete[] siteMap;
}

void inline SfluidRotationalObservable::display()
{

}

int inline SfluidRotationalObservable::read()
{
	return 0;
}

int inline SfluidRotationalObservable::write()
{
	return 0;
}

int inline SfluidRotationalObservable::readViaIndex(int idx)
{
	return 0;
}

int inline SfluidRotationalObservable::writeViaIndex(int idx)
{
	double iZ = 1.0/Zcount;
	double multC = 0.5/state->beta*iZ;
	Zcount = 0;

	stringstream ss;
	ss<<idx;
	string tfname = this->fileName + ss.str();

	ofstream wfp(tfname);
	wfp.precision(FIELDPRECISION);
	wfp.width(FIELDWIDTH);
	wfp.setf(FIELDFORMAT);

	double rhoS_x,rhoS_y,rhoS_z;
	double rhoA_x,rhoA_y,rhoA_z;
	rhoS_x = rhoS_y = rhoS_z = 0.0;
	rhoA_x = rhoA_y = rhoA_z = 0.0;

	for(int i=0;i<size;i++)
	{
		//Find coordinates relative to center of trap
		int z = siteMap[3*i];
		int y = siteMap[3*i+1];
		int x = siteMap[3*i+2];

		int rz2 = y*y+x*x;
		int ry2 = x*x+z*z;
		int rx2 = y*y+z*z;

		rhoS_z += rho_s[6*i];
		rhoA_z += rho_s[6*i+1];

		rhoS_y += rho_s[6*i+2];
		rhoA_y += rho_s[6*i+3];

		rhoS_x += rho_s[6*i+4];
		rhoA_x += rho_s[6*i+5];

		rho_s[6*i] = (rz2==0) ? 0.0: rho_s[6*i]/rz2*multC;
		rho_s[6*i + 1] = (rz2==0) ? 0.0: rho_s[6*i+1]/rz2*multC;

		rho_s[6*i + 2] = (ry2==0) ? 0.0: rho_s[6*i+2]/ry2*multC;
		rho_s[6*i + 3] = (ry2==0) ? 0.0: rho_s[6*i+3]/ry2*multC;

		rho_s[6*i + 4] = (rx2==0) ? 0.0: rho_s[6*i+4]/rx2*multC;
		rho_s[6*i + 5] = (rx2==0) ? 0.0: rho_s[6*i+5]/rx2*multC;

		wfp<<rho_s[6*i]<<" "<<rho_s[6*i+1]<<" "<<rho_s[6*i+2]<<" "<<rho_s[6*i+3]<<" "<<rho_s[6*i+4]<<" "<<rho_s[6*i+5]<<endl;

		rho_s[6*i] = rho_s[6*i+1] = rho_s[6*i+2] = rho_s[6*i+3] = rho_s[6*i+4] = rho_s[6*i+5] = 0.0;
	}
	wfp.close();


	//Write superfluid density
	wfp.open(this->fileName,std::ofstream::app);
	wfp<<rhoS_z*multC<<" "<<rhoA_z*multC<<" "<<rhoS_y*multC<<" "<<rhoA_y*multC<<" "<<rhoS_x*multC<<" "<<rhoA_x*multC<<" "<<endl;
	wfp.close();
}

void inline SfluidRotationalObservable::measure()
{

	for(int i=0;i<state->offdiagOperatorLocs.size();i++)
	{
		unsigned int opstate = state->opString[state->offdiagOperatorLocs[i]];
		int posloc = opstate%hamiltonian->BT - 1;
		//0 means diagonal 1,2 are off-diagonal entries; but here it is always 1 or 2 since no diagonal
		//Propagate  opType = 1 => a(^+)a, opType = 2 => aa(^+).
		//Acts on |s_k s_kp>
		//Corresponding posloc = 0 and 1

		unsigned int bond = opstate/hamiltonian->BT;
		int tt = (posloc == 1) ? 1 : -1;

		long boffset = bond*hamiltonian->bondMapOffset;
		int s_kp = 3*hamiltonian->bondMap[boffset + 1-posloc];
		int s_k = 3*hamiltonian->bondMap[boffset + posloc];

		//Find real coordinates relative to specified axis
		int zk = siteMap[s_k];
		int yk = siteMap[s_k+1];
		int xk = siteMap[s_k+2];

		int zkp = siteMap[s_kp];
		int ykp = siteMap[s_kp+1];
		int xkp = siteMap[s_kp+2];

		//Z-axis
		double lval = (yk*xkp - ykp*xk);
		Atotal_r[s_k] += lval;
		Atotal[0] += lval;

		//Y-axis
		lval = (zk*xkp - zkp*xk);
		Atotal_r[s_k+1] += lval;
		Atotal[1] += lval;

		//X-axis
		lval = (yk*zkp - ykp*zk);
		Atotal_r[s_k+2] += lval;
		Atotal[2] += lval;

	}

	for(int ss=0;ss<size;ss++)
	{
		double val0 = Atotal_r[3*ss];
		double val1 = Atotal_r[3*ss+1];
		double val2 = Atotal_r[3*ss+2];

		rho_s[6*ss] += Atotal[0]*val0;
		rho_s[6*ss+1] += val0*val0;
		rho_s[6*ss+2] += Atotal[1]*val1;
		rho_s[6*ss+3] += val1*val1;
		rho_s[6*ss+4] += Atotal[2]*val2;
		rho_s[6*ss+5] += val2*val2;

		Atotal_r[3*ss] = Atotal_r[3*ss+1] = Atotal_r[3*ss+2] = 0.0;
	}

	Atotal[0] = Atotal[1] = Atotal[2] = 0.0;
	Zcount++;
}

void inline SfluidRotationalObservable::clear()
{
	for(int i=0;i<3*size;i++)
		Atotal_r[i] = 0.0;

	for(int i=0;i<6*size;i++)
			rho_s[i] = 0.0;

	Zcount = 0;
	Atotal[0] = Atotal[1] = Atotal[2] = 0.0;
}

int inline SfluidRotationalObservable::checkPointWrite()
{
	ofstream wif(chkFileName);
	wif.precision(FIELDPRECISION);
	wif.width(FIELDWIDTH);
	wif.setf(FIELDFORMAT);

	if(wif)
	{
		wif<<Zcount<<endl;
		for(int i=0;i<size;i++)
			wif<<rho_s[6*i]<<" "<<rho_s[6*i+1]<<" "<<rho_s[6*i+2]<<" "<<rho_s[6*i+3]<<" "<<rho_s[6*i+4]<<" "<<rho_s[6*i+5]<<endl;
	}
	wif.close();
	return 0;
}

int inline SfluidRotationalObservable::checkPointRead()
{
	ifstream rif(chkFileName);
	if(rif)
	{
		rif>>Zcount;
		for(int i=0;i<size;i++)
			rif>>rho_s[6*i]>>rho_s[6*i+1]>>rho_s[6*i+2]>>rho_s[6*i+3]>>rho_s[6*i+4]>>rho_s[6*i+5];
	}
	rif.close();
	return 0;
}


} /* namespace measures */
