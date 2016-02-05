/*
 * GreenObservable.cpp
 *
 *  Created on: May 3, 2014
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

using namespace std;
using namespace __gnu_cxx;
using namespace measures;

GreenObservable::GreenObservable(string s, core::StateVariable* lstate,
core::Hamiltonian* lham, FILE* llog, Flags& f, GijI& _loopVar) : Measurable(s,lstate,lham,llog), lpvar(_loopVar)
{
	this->flag = f;
	Zcount = 0;
	this->sigWidth = flag.sigma;
	numTvalues = 0;
}

GreenObservable::~GreenObservable()
{
	tcollection.clear();
	numTvalues = 0;
	Zcount = 0;

	orderedMap2d::iterator ita = tvalues.begin();
	for(;ita!=tvalues.end();++ita)
			ita->second.clear();
	tvalues.clear();

	clear();
}

void GreenObservable::createTvalues()
{
	L = state->opCutOff;

	if(flag.fixedGreenTime)
	{
		tvalues[0][0] = 1.0;
		numTvalues++;
		tcollection.push_back(0.0);

		double tm = flag.greendt;
		do
		{
			tcollection.push_back(tm);
			double sig = sigWidth*sqrt(L*tm*(1.0-tm));
			int cl = int(L*tm);
			int sigw = int(sig);
			int lb = floor(cl-sig); lb = (lb<0) ? 0: lb;
			int ub = floor(cl+sig); ub = (ub>=L) ? L-1:ub;

			for(int i=lb;i<=ub;i++)
				tvalues[numTvalues][i] = gsl_cdf_binomial_P(i,tm,L);

			numTvalues++;
			tm+=flag.greendt;
		}while(tm<1.0 && abs(tm-1.0)>1.0e-6);
	}
	else
	{
		double s2 = sigWidth*sigWidth;

		//Generate t points
		tvalues[0][0] = 1.0; // t = 0.0, \delta_l = 0
		numTvalues++;
		tcollection.push_back(0.0);

		double tm = 0.0;
		while(tm<1.0)
		{
			double c = L*tm + sigWidth*sqrt(L*tm*(1.0-tm));

			double a = 2*L*(L+s2);
			double b = L*(2*c+s2);
			double tmn = (b + sqrt(b*b-2*a*c*c))/a;

			if(abs(tmn-tm)<1.0e-9) //ZEROTOL
				break;

			tm = tmn;
			tcollection.push_back(tm);

			double sig = sigWidth*sqrt(L*tm*(1.0-tm));
			double cl = L*tm;
			int lb = floor(cl-sig)+1;
			int ub = floor(cl+sig);

			lb = (lb<0) ? 0: lb;
			ub = (ub>=L) ? L-1:ub;
			for(int i=lb;i<=ub;i++)
				tvalues[numTvalues][i] = gsl_cdf_binomial_P(i,tm,L);

			numTvalues++;
		}
	}

	//Check if output specified
	if(flag.outputTimeSeries)
	{

		ofstream of(flag.timeSeriesFile);

		of<<"Maximum level: "<<L<<endl;

		for(int i=0;i<tcollection.size();i++)
		{
			sortedRow::iterator itb = tvalues[i].begin();

			sortedRow::iterator ite = (--tvalues[i].end());

			of<<i<<"\t"<<tcollection[i]<<"\t\t"<<tvalues[i].size()
					<<"\t\t"<<itb->first<<"\t\t"<<itb->second<<"\t\t"
					<<ite->first<<"\t\t"<<ite->second<<endl;
		}
		of.close();
	}
}

void GreenObservable::measure()
{
	Zcount++;


	GijI::iterator it = lpvar.begin();
	for(;it!=lpvar.end();++it)
	{
		int sitei = it->first;
		SiteToInterval::iterator it2 = it->second.begin();
		for(;it2!=it->second.end();++it2)
		{
			int sitej = it2->first;

			for(int tidx=0;tidx<numTvalues;tidx++)
			{
				//Get lower and upper bounds
				int flb = tvalues[tidx].begin()->first;
				int fub = (--tvalues[tidx].end())->first;

				double sumc  = 0.0;
				Interval::iterator it3 = it2->second.begin();
				for(;it3!=it2->second.end();++it3)
				{
					int delta_l = it3->first;
					int delta_u = it3->second.endLevel;
					float val = it3->second.value;

					//Must handle 0 separately
					if(delta_l==0 && delta_u == 0)
					{
						if(tvalues[tidx].find(0)!=tvalues[tidx].end())
							sumc += val*tvalues[tidx][0];
						continue;
					}

					//Now since 0 has been handled we can take normal boundary interval spec
					delta_l--; //Recall I = ]delta_l1,delta_l2]

					if(fub<=delta_l || delta_u<flb)
						continue; //Convolution will be 0
					else if(delta_u==flb)
					{
						sumc +=	val*tvalues[tidx][flb];
						continue;
					}

					int lptr,uptr; //Get relevant bounds
					if(delta_l<=flb && fub<=delta_u)
					{
						lptr = flb;
						uptr = fub;
					}
					else if(flb<delta_l && delta_u<fub)
					{
						lptr = delta_l;
						uptr = delta_u;
					}
					else if(flb<delta_l && fub>delta_l)
					{
						lptr = delta_l;
						uptr = fub;
					}
					else if(flb<delta_u && fub>delta_u)
					{
						lptr=flb;
						uptr=delta_u;
					}

					sumc +=	val*(tvalues[tidx][uptr] - tvalues[tidx][lptr]);
				}

				if(sumc>0.0)
					green[sitei][sitej][tidx] = green[sitei][sitej][tidx] + sumc;
			}
			it2->second.clear(); //Once time stream is done clear out intervals
		}
		it->second.clear();
	}
	lpvar.clear();
}

void GreenObservable::clear()
{
	map3d::iterator it = green.begin();
	for(;it!=green.end();++it)
	{
		map2d::iterator it2 = it->second.begin();
		for(;it2!=it->second.end();++it2)
			it2->second.clear();
		it->second.clear();
	}
	green.clear();
}

void GreenObservable::display()
{

}

int GreenObservable::write()
{

}

int GreenObservable::read()
{

}

int GreenObservable::writeViaIndex(int idx)
{
	double iZ = 1.0/Zcount/state->lpContrIter;
	Zcount = 0;

	stringstream ss;
	ss<<idx;
	string tfname = this->fileName + ss.str();

	ofstream wfp(tfname);
	wfp.precision(FIELDPRECISION);
	wfp.width(FIELDWIDTH);
	wfp.setf(FIELDFORMAT);

	map3d::iterator it = green.begin();
	for(;it!=green.end();++it)
	{
		int sitei = it->first;
		map2d::iterator it2 = it->second.begin();
		for(;it2!=it->second.end();++it2)
		{
			int sitej = it2->first;

			Row::iterator it3 = it2->second.begin();
			for(;it3!=it2->second.end();++it3)
			{
				int tval = it3->first;
				double value = it3->second*iZ;

				wfp<<sitei<<" "<<sitej<<" "<<tcollection[tval]<<" "<<value<<endl;
			}
		}
	}

	wfp.close();
	clear();
}

int GreenObservable::readViaIndex(int idx)
{

}

int GreenObservable::checkPointWrite()
{
	ofstream wfp(chkFileName);
	wfp.precision(FIELDPRECISION);
	wfp.width(FIELDWIDTH);
	wfp.setf(FIELDFORMAT);

	map3d::iterator it = green.begin();
	for(;it!=green.end();++it)
	{
		int sitei = it->first;
		map2d::iterator it2 = it->second.begin();
		for(;it2!=it->second.end();++it2)
		{
			int sitej = it2->first;

			Row::iterator it3 = it2->second.begin();
			for(;it3!=it2->second.end();++it3)
			{
				int tval = it3->first;
				double value = it3->second;

				wfp<<sitei<<" "<<sitej<<" "<<tval<<" "<<value<<endl;
			}
		}
	}

	wfp.close();
}

int GreenObservable::checkPointRead()
{
	int i,j,t;
	double val;

	ifstream rfp(chkFileName);
	while(rfp)
	{
		rfp>>i>>j>>t>>val;
		green[i][j][t] = val;
	}
	rfp.close();
}


