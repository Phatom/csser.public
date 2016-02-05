/*
 * TransitionMap.cpp
 *
 *  Created on: Apr 21, 2014
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
#include "SSTransitionMap.h"

using namespace std;
using namespace __gnu_cxx;
using namespace core;

SSTransitionMap::~SSTransitionMap()
{
	if(this->mapsInitialized)
	{
		delete[] opmap;
		delete[] statemap;
		delete[] vtxTypeMap;
		delete[] newVtxMap;
		mapsInitialized = false;
	}

	if(this->probsCalculated)
	{
		delete[] bondToProb;
		probsCalculated = false;
	}
}

void SSTransitionMap::initializeMaps()
{
	//Already initialized so do nothing
	if(mapsInitialized)
		return;

	VM = (hamiltonian->LAM+1);
	VM2 = VM*VM;
	VM3 = (hamiltonian->LAM+1)*VM2;
	VM4 = (hamiltonian->LAM+1)*VM3;

	long totalAS = (hamiltonian->LAM-1)*(3*hamiltonian->LAM+1) + 2*(2*hamiltonian->LAM+1);
	opmap = new char[4*4*2];
	statemap = new char*[totalAS]; //Given vertex type what are the states of the four legs
	vtxTypeMap = new long[VM4]; //Given  state of 4 legs what is vertex type

	//For operator types
	for(int op=0;op<2;op++)
		for(int enll=0;enll<4;enll++)
			for(int exll=0;exll<4;exll++)
			{
				if(enll/2 > 0 && exll/2 > 0 || enll/2 < 1 && exll/2 < 1)
					opmap[8*enll + 2*exll + op] = 1 - op;
				else
					opmap[8*enll + 2*exll + op] = op;
			}

	vtxCount = 0;
	for(int i=0,ii=0;i<VM;i++,ii+=VM3)
		for(int j=0,jj=0;j<VM;j++,jj+=VM2)
			for(int k=0,kk=0;k<VM;k++,kk+=VM)
				for(int l=0;l<VM;l++)
				{
					if(i==k && j==l || i==k+1 && j==l-1 || i==k-1 && j==l+1)
					{
							vtxTypeMap[ii+jj+kk+l] = vtxCount;

							statemap[vtxCount] = new char[4];
							statemap[vtxCount][0] = i;
							statemap[vtxCount][1] = j;
							statemap[vtxCount][2] = k;
							statemap[vtxCount][3] = l;

							vtxCount++;


#ifdef ERRORCHECK
							if(p>totalAS)
								break;
#endif
					}
					else
					{
							vtxTypeMap[ii+jj+kk+l] = -1;
					}
				}

	fprintf(log,"Allocated size for acceptable states: %d\n",totalAS);
	fprintf(log,"Actual no. of acceptable states: %d\n",vtxCount);
	fflush(log);

	//New Vertex Type Map
	vtxOffset = vtxCount*32;
	newVtxMap = new long[vtxOffset]; //Given vertex type, entrance leg no. and exit leg no., and operator what is new vertex

	int newStates[4];
	const char* currStates;
	char newOp;

	for(unsigned long i=0, ii=0;i<vtxCount;i++,ii+=32) //Vertex
		for(int j=0, jj=0;j<4;j++,jj+=8) //Entrance
			for(int k=0, kk=0;k<4;k++,kk+=2) //Exit
				for(int o=0;o<2;o++) //Operator
				{
					//if(statemap[i][0] == -1)
					//	newVtxMap[i][j][k][o]=-1;//Vertex not possible
					if(j==k)
					{       //Bounces are NOT GUARANTEED TO WORK ALL THE TIME
						currStates = statemap[i];
						char nsj1 = currStates[j];
						char nsj2 = currStates[k];

						newOp = opmap[8*j + 2*k + o];

						nsj1 += ((o == 0) ? 1:-1);
						//nsj2 += ((newOp == 0) ? 1:-1);

						if(nsj1 >= (hamiltonian->LAM+1)  || nsj1<0)
							newVtxMap[ii+jj+kk+o] = -1;
						else
							newVtxMap[ii+jj+kk+o] = i;
					}
					else
					{
						currStates = statemap[i]; //Get current states
						newStates[0] = currStates[0];
						newStates[1] = currStates[1];
						newStates[2] = currStates[2];
						newStates[3] = currStates[3];

						newOp = opmap[8*j + 2*k + o]; //Get new operator

						//0 = a+  1 = a
						newStates[j] += ((o == 0) ? 1: -1);
						newStates[k] += ((newOp == 0) ? 1: -1);

						if(newStates[j] >= VM || newStates[k] >= VM || newStates[j]<0 || newStates[k] < 0) //Either full or empty
							newVtxMap[ii+jj+kk+o] = -1; //vertex not possible
						else
							newVtxMap[ii+jj+kk+o]= vtxTypeMap[ VM3*newStates[0] + VM2*newStates[1] + VM*newStates[2] + newStates[3] ];
					}
				}

	this->mapsInitialized = true;
}

double SSTransitionMap::calcExP(int bond, int si, int sj, long vtxState, char ladderOp, char entranceLeg, double* tmatrix)
{

	long vVCount = 0;

	long smap[4][2] = { {-1, -1}, {-1, -1}, {-1, -1}, {-1, -1} } ;
	char vmap[4][4];
	double a[4+6+2][16+1];
	double weights[4][4];

	long i = vtxState;
	char enl = entranceLeg;
	char o = ladderOp;

	char j,k;
	for(j=0;j<12;j++)
		for(k=0;k<17;k++)
			a[j][k] = 0.0;

	for(j=0;j<4;j++)
		for(k=0;k<2;k++)
			smap[j][k]=-1;

	if(newVtxMap[32*i + 8*enl + 2*enl + o] != -1)
	{
		smap[enl][0] = i;
		smap[enl][1] = o;
	}
	else
	{
		smap[enl][0] = -1;
		smap[enl][1] = -1;
	}


	for(k=0;k<4;k++)
	{
		if(smap[enl][0] != -1 && enl!=k)
		{
#ifdef ERRORCHECK
			if(smap[k][0] != -1 && (smap[k][0]!=newVtxMap[ 32*smap[enl][0] + 8*enl + 2*k + smap[enl][1] ] || smap[k][1] != 1 - opmap[8*enl + 2*k + smap[enl][1]]))
			{
				return 1.0;
				fprintf(stdo,"Map Issues!\n");
			}
#endif
			//Update next/previous
			smap[k][0] = newVtxMap[32*smap[enl][0]+ 8*enl + 2*k + smap[enl][1]];
			smap[k][1] = 1 - opmap[8*enl + 2*k + smap[enl][1]]; //Operator is conjugate
		}
	}

	//Calculate Weights
	for(int j=0;j<4;j++) //Current entrance leg was this
		for(int k=0;k<4;k++) //Possible Exit Leg
		{
			if(smap[j][0] == -1)
			{
				weights[j][k] = 0.0;
				continue;
			}

			const char* newStates = statemap[smap[j][0]]; //Find new leg states
			char entlegState = newStates[j];

			//Remember operator 0 = a+ 1 = a
			//unsigned int newOpConj = 1 - newOp; //Find conjugate of new operator
			char currentOperator = smap[j][1];
			int sgn = (currentOperator == 0) ? 1: 0;

			double factor = hamiltonian->operateWith(currentOperator,entlegState); //sqrt(entlegState + sgn);
			//double factor = 1.0;
			char idiff = newStates[0]-newStates[2];
			char jdiff = newStates[1]-newStates[3];

			if(idiff==0 && jdiff==0) //diagonal
			{
				weights[j][k] = factor*hamiltonian->getDiagonalTerm(si,sj,newStates[0],newStates[1]);
				//(global.C - 0.5/global.Dimension*(E_i*newStates[0] + E_j*newStates[1] - global.mu*(newStates[0] + newStates[1]) + (U_i*newStates[0]*(newStates[0]-1) + U_j*newStates[1]*(newStates[1]-1))) + global.V*newStates[0]*newStates[1]);
			}
			else //off-diagonal
			{
				weights[j][k] = factor*hamiltonian->getOffdiagonalTerm(bond,currentOperator,newStates);
			}

		}

	//Setup variables
	char p = 1;
	for(j=0;j<4;j++)
	{
		for(k=0;k<4;k++)
		{
			if(smap[j][0]!=-1 && smap[k][0]!=-1) //Transition allowed
				vmap[j][k] = p++;
			else
				vmap[j][k] = -1;
//			cout<<(int) vmap[j][k]<<" ";
		}
//		cout<<endl;
	}

	int N = p-1; //Actual number of variables

	a[0][0] = 0.0; //Minimization requirement
	p = 1;
	//Sum = 1 constraint
	for(j=0;j<4;j++)
	{
		a[p][0] = 1.0;
		char flag = 0;
		for(k=0;k<4;k++)
		{
			//If transition is allowed then variable will exist.
			if(vmap[j][k]!=-1)
			{
				flag = 1;
				a[p][vmap[j][k]] = -1.0;
			}
		}

		//Minimization only makes sense if actual transitions are allowed
		if(flag)
		{
			a[0][vmap[j][j]] = -1.0;
			p++;
		}
	}

	//Actual number of equations;
	int M = p - 1;

	//cout<<"Detailed Balance"<<endl;
	//Detailed Balance constraint
	//Note that we add the k->k transition as well.
	char s=1;
	for(j=0;j<3;j++)
	{
		for(k=s;k<4;k++)
		{
			if(vmap[j][k]!=-1)
			{
				a[p][0] = 0.0;
				a[p][vmap[j][k]] = weights[j][k];
				a[p][vmap[k][j]] = -weights[k][j];
				p++;
			}
		}
		s++;
	}

	M = p - 1; //Actual number of equations

	Mat_DP mata(M+2,N+1);
	for(j=0;j<M+2;j++)
	{
		for(k=0;k<N+1;k++)
		{
			mata[j][k] = (double) a[j][k];
		}
	}

	Vec_INT izrov(N), iposv(M);

/*
				//OUTPUT
				cout<<"Pre output"<<endl;
				for(j=0;j<M+2;j++)
				{
					for(k=0;k<N+1;k++)
					{
						printf("% 3.2lf ",mata[j][k]);
					}
					cout<<endl;
				}
				cout<<"$$$$"<<endl;
*/

	int scase;
	if(M>0 && N>0)
		NR::simplx(mata,0,0,M,scase,izrov,iposv);

	if(scase == -7)
	{
		fprintf(log,"SOLUTION DID NOT CONVERGE WITHIN SPECIFIED LIMIT!\n");
		return 1.0;
	}

	int n=0;
	for(j=0;j<4;j++)
	{
//		if(smap[j][0]!=-1)
//			fprintf(stdo,"% 5ld[% d % d % d % d] % 5ld ",smap[j][0],statemap[smap[j][0]][0],statemap[smap[j][0]][1],statemap[smap[j][0]][2],statemap[smap[j][0]][3],smap[j][1]);
//		else
//			fprintf(stdo,"% 5ld[% d % d % d % d] % 5ld ",smap[j][0],-1,-1,-1,-1,smap[j][1]);
		for(k=0;k<4;k++)
		{
//
//			if(vmap[j][k]==-1)
//				fprintf(stdo,"% 3.2lf {HUH?} ",0.0);


			char flag=0;
			for(p=0;p<M;p++)
			{
				if(iposv[p]+1==vmap[j][k])
				{
					if(mata[p+1][0]<0.0)
						mata[p+1][0] = 0.0;

					tmatrix[32*smap[j][0] + 8*j + smap[j][1] + 2*k] = mata[p+1][0];
					flag=1;
//					fprintf(stdo,"% 3.2lf ",eProbs[4*j+k]);
					break;
				}
			}

//			if(!flag && vmap[j][k]!=-1)
//				fprintf(stdo,"% 3.2lf ",0.0);
		}
//		cout<<endl;
	}
//	cout<<"=================================="<<endl<<endl;

	return -mata[0][0];
}

void SSTransitionMap::calcExitProbs(string base)
{
	if(this->probsCalculated)
	{
		//If already calculated and being asked to do so again
		//Means some parameters must have been changed
		//So cleanup and then re-calculate

		delete[] bondToProb;
		probsCalculated = false;
	}

	fprintf(log,"==========================================================================\n");
	fprintf(log,"Transition Probability Stats\n");
	fprintf(log,"==========================================================================\n");

	//Store bond to probabilities
	//This is memory intensive
	this->tsize = (hamiltonian->symmetricSetSize)*this->vtxOffset;
	bondToProb = new double[tsize];
	//BELOW IS VERY IMPORTANT
	for(int i=0;i<tsize;i++)
		bondToProb[i] = 0.0;

	if(readBondProbs(base) == SUCCESS)
	{
		fprintf(log,"Transition probabilities loaded from file.\n\n");
		probsCalculated = true;
		return;
	}

	///////////////////////////////////////////////////////////////////
	// CALCULATE TRANSITION PROBABILITIES
	///////////////////////////////////////////////////////////////////

	int si,sj; //Sites
	int x,xp,y,yp,z,zp;

	///////////////////////////////////////////////////////////////////
	double maxbp = -1.0;
	int maxsi,maxsj;
	int maxState, maxop, maxenl;
	///////////////////////////////////////////////////////////////////

	long p;
	#pragma omp parallel for default(none) shared(maxbp,maxsi,maxsj,maxState,maxop,maxenl) private(p)
	for(p=0;p<hamiltonian->symmetricSetSize;p++)
	{
		int loc = hamiltonian->symmSetOffset*p;

		int si = hamiltonian->symmetricSet[loc];
		int sj = hamiltonian->symmetricSet[loc + 1];
		int bond = hamiltonian->pairSitesToBond[pair<int,int>(si,sj)];

		double lmaxbp;
		int lmaxsi,lmaxsj,lmaxState,lmaxop,lmaxenl;

		for(int vtx=0;vtx<vtxCount;vtx++)
			for(char op=0;op<2;op++)
				for(char enl=0;enl<4;enl++)
				{
					double bounceProb = calcExP(bond,si,sj,vtx,op,enl,bondToProb+vtxOffset*p);
						
					if(bounceProb >= lmaxbp)
					{
						lmaxbp = bounceProb;
						lmaxsi = si;
						lmaxsj = sj;
						lmaxState = vtx;
						lmaxop = op;
						lmaxenl = enl;
					}
				}

		#pragma omp critical
		{
			if(lmaxbp >= maxbp)
			{

				maxbp = lmaxbp;
				maxsi = lmaxsi;
				maxsj = lmaxsj;
				maxState = lmaxState;
				maxop = lmaxop;
				maxenl = lmaxenl;
			}
		}
	}

	fprintf(log,"Max. Bounce: %6.4e\n",maxbp);
	fprintf(log,"State: %d %d %d %d\n",statemap[maxState][0],statemap[maxState][1],statemap[maxState][2],statemap[maxState][3]);
	fprintf(log,"Sites: %d,%d\n",maxsi,maxsj);
	fprintf(log,"Operator: %d\n",maxop);
	fprintf(log,"Entrance Leg: %d\n",maxenl);
	fflush(log);

	//Now calculate cumulatives
#pragma omp parallel for default(none) private(p)
	for(p=0;p<hamiltonian->symmetricSetSize;p++)
	{
		int loc = hamiltonian->symmSetOffset*p;
		int si = hamiltonian->symmetricSet[loc];
		int sj = hamiltonian->symmetricSet[loc + 1];

		for(int vtx=0;vtx<vtxCount;vtx++)
		{
			for(char op=0;op<2;op++)
			{
				for(char enl=0;enl<4;enl++)
				{
					for(char exl=1;exl<4;exl++)
					{
						unsigned int pr = vtxOffset*p + vtx*32 + 8*(enl) + 2*exl + op;
						unsigned int prm = vtxOffset*p + vtx*32 + 8*(enl) + 2*(exl-1) + op;

						bondToProb[pr] += bondToProb[prm];
					}

					//Round-off issues are important
					unsigned int pr = vtxOffset*p + vtx*32 + 8*(enl) + 2*3 + op;
					if(fabs(bondToProb[pr])-1.0>1.0e-10)
					{
						fprintf(log,"Transition setup issues!\n");
						for(int ll=0;ll<4;ll++)
							fprintf(log,"Vtx: %d(%d %d %d %d)  Enl: %d  Operator: %d  Exit leg: %d  Prob:  %9.6le\n",vtx,
								statemap[vtx][0],
								statemap[vtx][1],
								statemap[vtx][2],
								statemap[vtx][3],
								enl,op,ll,bondToProb[pr+ll*2 - 6]);
						fprintf(log,"------------------------------------------------------------------------------\n");
					}
					
					if(fabs(bondToProb[pr])>1.0e-10) //i.e. not zero
						bondToProb[pr] = 1.0;
				}
			}
		}
	}

	///////////////////////////////////////////////////////////////////
	// WRITE PROBABILITIES
	///////////////////////////////////////////////////////////////////
	writeBondProbs(base);
	probsCalculated = true;
	fprintf(log,"Transition probabilities written to file.\n\n");
	///////////////////////////////////////////////////////////////////

#if 0
	/*
	FILE* fp = fopen(DEBUGOUT,"w");

	for(it = rpairMap.begin();it!=rpairMap.end();++it)
	{
		fprintf(fp,"(%d,%d) %d\n",it->first.first,it->first.second,it->second);
	}
	fprintf(fp,"\n\n");

	for(long di=0;di<Nb;di++)
		fprintf(fp,"bond[%ld] = %d\n",di,bondToRpair[di]);
		*/

	for(p=0;p<hamiltonian->symmetricSetSize;p++)
	{
		int loc = hamiltonian->symmSetOffset*p;
		int si = hamiltonian->symmetricSet[loc];
		int sj = hamiltonian->symmetricSet[loc + 1];
		fprintf(log,"%d = %d %d\n",loc,si,sj);

		for(int vtx=0;vtx<vtxCount;vtx++)
		{
			for(char op=0;op<2;op++)
			{
				for(char enl=0;enl<4;enl++)
				{
					for(char exl=0;exl<4;exl++)
					{
						unsigned int pr = vtxOffset*p + vtx*32 + 8*(enl) + 2*exl + op;

						fprintf(log,"[%d %d %d %d] %d %d %d %10.6e\n",
								statemap[vtx][0],
								statemap[vtx][1],
								statemap[vtx][2],
								statemap[vtx][3],
								op,(int) enl, (int) exl,
								bondToProb[pr]);
					}
					fprintf(log,"\n");
				}
			}
		}


	}

#endif
}
