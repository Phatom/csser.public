/*
 * BasicSimEngine.cpp
 *
 *  Created on: Apr 21, 2014
 *      Author: ushnish
 *
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

using namespace std;
using namespace __gnu_cxx;
using namespace core;


long SPDMSimEngine::loopUpdate(long long& lup_er)
{
	//Not using vtxToBond but rather vtxToSG
	unsigned int* vtxToBond = new unsigned int[stateVariable->nOp];
	unsigned int* vtxToSG = new unsigned int[stateVariable->nOp];

	memcpy(state_t,stateVariable->state,hamiltonian->numSites*sizeof(char));

	for(int j=0;j<hamiltonian->numSites;j++)
		first[j] = last[j]  = -1;

	if(stateVariable->nOp>0)
	{
		unsigned int* optovtx = new unsigned int[stateVariable->opCutOff];
		unsigned int p = 0, b;
		unsigned int s0,s1;

		for(unsigned int j=0;j<stateVariable->opCutOff;j++)
		{
			if(stateVariable->opString[j]!=0)
			{
				unsigned int poffset = 4*p;

				b = stateVariable->opString[j]/hamiltonian->BT;
				unsigned int boffset = hamiltonian->bondMapOffset*b;

				s0 = hamiltonian->bondMap[boffset];
				s1 = hamiltonian->bondMap[boffset + 1];

				//Get vtx to bond map
				vtxToBond[p] = b;
				//Get vtx to symmetric group map
				vtxToSG[p] = hamiltonian->bondToSSet[b];

				optovtx[j] = p;

				if(last[s0]!=-1)
				{
					stateVariable->link[poffset] = last[s0];
					stateVariable->link[last[s0]] = poffset;
					last[s0] = poffset + 2;
				}
				else if(last[s0]==-1)
				{
					first[s0] = poffset;
					last[s0] = poffset + 2;
				}

				if(last[s1]!=-1)
				{
					stateVariable->link[poffset+1] = last[s1];
					stateVariable->link[last[s1]] = poffset + 1;
					last[s1] = poffset + 3;
				}
				else if(last[s1]==-1)
				{
					first[s1] = poffset + 1;
					last[s1] = poffset + 3;
				}

				char oldstate[2] = {state_t[s0],state_t[s1]};
				char optype = stateVariable->opString[j]%hamiltonian->BT;
				if(optype != 0) //Off-diagonal term
				{
					//Propagate  opType = 1 => a(^+)a, opType = 2 => aa(^+).
					char uvi = (optype == 1) ? 1: -1;
					state_t[hamiltonian->bondMap[boffset]] += uvi;
					state_t[hamiltonian->bondMap[boffset+1]] -= uvi;
				}

				stateVariable->vtx[p] = tmap->vtxTypeMap[tmap->VM3*oldstate[0] + tmap->VM2*oldstate[1] + tmap->VM*state_t[s0] + state_t[s1]];
				p++;
			}
		}

		for(unsigned int j=0;j<hamiltonian->numSites;j++)
		{
			if(last[j]!=-1)
			{
				stateVariable->link[last[j]] = first[j];
				stateVariable->link[first[j]] = last[j];
			}
		}


		for(long j=0;j<stateVariable->lpContrIter;j++)
		{
			int crossCheck = -1;
			long long localLoopLength = 0;
			unsigned int vtx0 = gsl_rng_uniform_int(rgenref,4*stateVariable->nOp);
			char enl_state = tmap->statemap[stateVariable->vtx[vtx0/4]][vtx0%4];

			char ladderOp;
			if(gsl_rng_uniform(rgenref)>=0.5)
				ladderOp = 1;
			else
				ladderOp = 0;

			if(enl_state==hamiltonian->LAM && ladderOp==0 || enl_state==0 && ladderOp==1)
				continue;

			unsigned int vtxC = vtx0;
			//char ladderOp0 = ladderOp;
			unsigned int ovtx0 = stateVariable->vtx[vtxC/4];
			unsigned int vtxp0 = vtx0/4; //Vertex number 0
			char enl0 = vtx0%4; //entrance leg 0

			char ovtx0state = tmap->statemap[stateVariable->vtx[vtxp0]][enl0]; //Number of particles at leg entrance
			char ovtx0measure = ovtx0state + 1 - ladderOp; //Crosses will measure this
			int finalSite = hamiltonian->bondMap[hamiltonian->bondMapOffset*vtxToBond[vtxp0]+enl0%2]; //Actual site information

			do{
				p = vtxC/4; //Vertex Number
				char enl = vtxC%4; //Entrance leg

				unsigned int probptr = stateVariable->vtx[p]*32 + 8*enl + ladderOp; //map to map minus bond location
				double* entry = tmap->bondToProb;

				//NOTE THAT WE ARE USING SYMMETRIC GROUP AND NOT vtxToBond
				unsigned int offset = tmap->vtxOffset*vtxToSG[p] + probptr; //Full map to probability table
				//Find exit leg

				double rnd = gsl_rng_uniform(rgenref);
				char exl = (rnd<entry[offset]) ? 0 : ((rnd<entry[offset+2]) ? 1 : (rnd<entry[offset+4] ? 2 : (rnd<entry[offset+6] ? 3: -1)));

#if 0
				/*
				printf("rnd: %lf exl: %d\n",rnd,int(exl));
				for(int ll=0;ll<4;ll++)
					printf("Bond: %d  Vtx: %d(%d %d %d %d)  Enl: %d  Operator: %d  Exit leg: %d  Prob:  %9.6le\n",
							vtxToSG[p],stateVariable->vtx[p],
							tmap->statemap[stateVariable->vtx[p]][0],
							tmap->statemap[stateVariable->vtx[p]][1],
							tmap->statemap[stateVariable->vtx[p]][2],
							tmap->statemap[stateVariable->vtx[p]][3],
							enl,ladderOp,ll,entry[offset+ll*2]);
				fflush(0);
				*/

				if(exl==-1)
				{
					cout<<"Error! No match!"<<endl;
					cout<<"M: "<<stateVariable->opCutOff<<endl;

					int oldp = vtx0/4;
					int oldenl = vtx0%4;
					cout<<"Inital entry point stats: "<<vtx0<<" "<<oldp<<" "<<(int) enl_state<<" "<<ovtx0<<endl;
					unsigned int pr = tmap->vtxOffset*vtxToSG[oldp] + ovtx0*32 + 8*(oldenl) + ladderOp0;

					for(int ll=0;ll<4;ll++)
						printf("Bond: %d  Vtx: %d(%d %d %d %d)  Enl: %d  Operator: %d  Exit leg: %d  Prob:  %9.6le\n",
								vtxToSG[oldp],ovtx0,
								tmap->statemap[ovtx0][0],
								tmap->statemap[ovtx0][1],
								tmap->statemap[ovtx0][2],
								tmap->statemap[ovtx0][3],
								oldenl,ladderOp0,ll,entry[pr+ll*2]);

					cout<<"Current Entry point: "<<p<<endl;
					for(int ll=0;ll<4;ll++)
						printf("Bond: %d  Vtx: %d(%d %d %d %d)  Enl: %d  Operator: %d  Exit leg: %d  Prob:  %9.6le\n",
								vtxToSG[p],stateVariable->vtx[p],
								tmap->statemap[stateVariable->vtx[p]][0],
								tmap->statemap[stateVariable->vtx[p]][1],
								tmap->statemap[stateVariable->vtx[p]][2],
								tmap->statemap[stateVariable->vtx[p]][3],
								enl,ladderOp,ll,entry[offset+ll*2]);
					fflush(0);
					return EXITCODE;
				}
#endif

				///////////////////////////////////////////////////////////
				if(enl!=exl)
					localLoopLength++; //Bounces do not add to loop length

				probptr += 2*exl;
				int nv = tmap->newVtxMap[probptr]; //Get new vertex type //REMEMBER THIS HAPPENS BECAUSE OF FLIPPING

				stateVariable->vtx[p] = nv; //Update vertex type
				ladderOp = tmap->opmap[8*enl + 2*exl + ladderOp]; //Update ladder operator

				vtxC = 4*p + exl;
				if(vtxC==vtx0)
					crossCheck = 6;
				else if(vtxC!=vtx0) //If it's reached then no vertex update of next vertex is required
				{
					unsigned int vtxN = stateVariable->link[vtxC]; //Next vertex will be this
					long vtxpN = vtxN/4; //next vertex number
					char enlN = vtxN%4; //next entrance leg

					if(localLoopLength > maxLoopLength)
					{
						delete[] vtxToSG;
						delete[] optovtx;
						//All added values needs to be subtracted.
						return -LOOPLENGTHEXCEEDED;
					}

					//GREEN'S FUNCTION CALCULATIONS GO HERE
					//check for cross
					crossCheck = 0;
					if((p == vtxp0 && (exl>=2 && enl0>=2 || exl<=1 && enl0<=1)) || (vtxpN == vtxp0 && (enlN>=2 && enl0>=2 || enlN<=1 && enl0<=1))) //If new or current vertex is same as original vertex check to see that legs are at same level or not. If yes then cross has occurred
						crossCheck = 1;
					else if( (p < vtxp0 && vtxpN > vtxp0) || (p > vtxp0 && vtxpN < vtxp0)) //Border cross
					{
						if(vtxC < vtxN && exl>=2 && enlN<=1)
							crossCheck = 2;
						else if(vtxC > vtxN && exl<=1 && enlN>=2)
							crossCheck = 3;
					}
					else if( (p > vtxp0 && vtxpN > vtxp0) || (p < vtxp0 && vtxpN < vtxp0))
					{
						if(vtxC < vtxN && exl<=1 && enlN>=2)
							crossCheck = 4;
						else if(vtxC > vtxN && exl>=2 && enlN<=1)
							crossCheck = 5;
					}

					vtxC = vtxN; //Next vertex will be this
				}

				if(crossCheck>0)
				{
					unsigned int b3 = hamiltonian->bondMapOffset*vtxToBond[p];
					int initialSite = hamiltonian->bondMap[b3+exl%2]; //Actual site information

					//char vtxstate = tmap->statemap[stateVariable->vtx[p]][exl];
					//double vtxmeasure = sqrt(vtxstate + 1 - ladderOp);

					double matrixval = (initialSite == finalSite) ? ovtx0state : ovtx0measure;
					if(matrixval > 0.0)
						spdm[initialSite][finalSite] = spdm[initialSite][finalSite] + matrixval;
				}

			}while(vtx0 != vtxC);

			//OTHER DIRECT MEASUREMENTS GO HERE
		}


		stateVariable->offdiagOperatorLocs.clear(); //Keep track of off-diagonal elements which may be used by observables
		for(long j=0;j<stateVariable->opCutOff;j++)
		{
			if(stateVariable->opString[j]!=0)
			{
				unsigned int b = stateVariable->opString[j]/hamiltonian->BT;
				//Propagate  opType = 1 => a(^+)a, opType = 2 => aa(^+).
				char* states = tmap->statemap[stateVariable->vtx[optovtx[j]]];

				//This is correct! THANK GOD!!
				//optype = 0 remains that way
				//ni-nf = 1 becomes optype = 2 (i.e. nf-ni = -1)
				//ni-nf = -1 becomes optype = 1 (i.e. nf-ni = 1)
				int optype = states[0]-states[2];
				optype += (optype >= 0) ? optype : 2;
				stateVariable->opString[j] = hamiltonian->BT*b + optype;//Update bond type

				if(optype!=0)
					stateVariable->offdiagOperatorLocs.push_back(j);
			}
		}

		delete[] optovtx;
	}

	for(int s=0;s<hamiltonian->numSites;s++)
	{
		if(first[s] == -1)
		{
			double rnd = gsl_rng_uniform(rgenref);
			if(rnd>0.5)
				stateVariable->state[s] = (char) gsl_rng_uniform_int(rgenref,hamiltonian->LAM+1); //Free world line update for sites that are not affected but operator string.
		}
		else
		{
			long p = first[s]/4; //Find vertex
			char enl = first[s]%4; //Find entrance leg

			stateVariable->state[s] = tmap->statemap[stateVariable->vtx[p]][enl]; //Find new state given entrance leg and vertex type.
		}
	}

	delete[] vtxToBond;
	delete[] vtxToSG;
	return 0;
}


