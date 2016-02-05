/*
 * Hamiltonian.cpp
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
#include "Hamiltonian.h"

using namespace std;
using namespace __gnu_cxx;
using namespace core;

void Hamiltonian::latticeCreate()
{
	int p=0;
	int ns = this->numSites;
	int bondnum = 1; //0 us special as in null

	if(this->bc == PBC)
	{
		if(this->dimension == 1)
		{
			for(int x0=0;x0<nspd;x0++)
			{
				int x1 = (1 + x0)%nspd;
				bondMap[this->bondMapOffset*(p+1)] = x0;
				bondMap[this->bondMapOffset*(p+1) + 1] = x1;
				bondMap[this->bondMapOffset*(p+1) + 2] = (x0 == nspd-1) ? 1 : 0;
				pairSitesToBond.insert(pair<pair<int,int>,int>(pair<int,int>(bondMap[this->bondMapOffset*(p+1)],bondMap[this->bondMapOffset*(p+1)+1]),p+1));
				p++;
			}
		}
		else if(this->dimension == 2)
		{
			for(int y0=0;y0<nspd;y0++)
				for(int x0=0;x0<nspd;x0++)
				{
					int x1 = (1 + x0)%nspd;
					bondMap[this->bondMapOffset*(p+1)] = p;
					bondMap[this->bondMapOffset*(p+1) + 1] = x1 + y0*nspd;
					bondMap[this->bondMapOffset*(p+1) + 2] = (x0 == nspd-1) ? 1 : 0;
					pairSitesToBond.insert(pair<pair<int,int>,int>(pair<int,int>(bondMap[this->bondMapOffset*(p+1)],bondMap[this->bondMapOffset*(p+1)+1]),p+1));

					int y1 = (1 + y0)%nspd;
					bondMap[this->bondMapOffset*(p + 1 + ns)] = p;
					bondMap[this->bondMapOffset*(p + 1 + ns) + 1] = x0 + (y1)*nspd;
					bondMap[this->bondMapOffset*(p + 1 + ns) + 2] = (y0 == nspd-1) ? 2 : 0;
					pairSitesToBond.insert(pair<pair<int,int>,int>(pair<int,int>(bondMap[this->bondMapOffset*(p+1+ns)],bondMap[this->bondMapOffset*(p+1+ns)+1]),p+1+ns));

					p++; //Go to next site //This can be replaced with x0 + y0*nspd btw
				}
		}
		else if(this->dimension == 3)
		{
			for(int z0=0;z0<nspd;z0++)
				for(int y0=0;y0<nspd;y0++)
					for(int x0=0;x0<nspd;x0++)
					{
						int x1 = (1+x0)%nspd;
						bondMap[this->bondMapOffset*(p+1)] = p;
						bondMap[this->bondMapOffset*(p+1) + 1] = x1 + (y0)*nspd + (z0)*nspd*nspd;
						bondMap[this->bondMapOffset*(p+1) + 2] = (x0 == nspd-1) ? 1 : 0;
						pairSitesToBond.insert(pair<pair<int,int>,int>(pair<int,int>(bondMap[this->bondMapOffset*(p+1)],bondMap[this->bondMapOffset*(p+1)+1]),p+1));
						bondnum++;

						int y1 = (1+y0)%nspd;
						bondMap[this->bondMapOffset*(p + 1 + ns)] = p;
						bondMap[this->bondMapOffset*(p + 1 + ns) + 1] = x0 + (y1)*nspd + (z0)*nspd*nspd;
						bondMap[this->bondMapOffset*(p + 1 + ns) + 2] = (y0 == nspd-1) ? 2 : 0;
						pairSitesToBond.insert(pair<pair<int,int>,int>(pair<int,int>(bondMap[this->bondMapOffset*(p+1+ns)],bondMap[this->bondMapOffset*(p+1+ns)+1]),p+1+ns));
						bondnum++;

						int z1 = (1+z0)%nspd;
						bondMap[this->bondMapOffset*(p + 1 + 2*ns)] = p;
						bondMap[this->bondMapOffset*(p + 1 + 2*ns) + 1] = x0 + (y0)*nspd + (z1)*nspd*nspd;
						bondMap[this->bondMapOffset*(p + 1 + 2*ns) + 2] = (z0 == nspd-1) ? 3 : 0;
						pairSitesToBond.insert(pair<pair<int,int>,int>(pair<int,int>(bondMap[this->bondMapOffset*(p+1+2*ns)],bondMap[this->bondMapOffset*(p+1+2*ns)+1]),(p+1+2*ns)));
						bondnum++;

						p++;
					}
		}
	}
	else if(this->bc == HW)
	{
		if(this->dimension == 1)
		{
			for(int x0=0;x0<nspd;x0++)
			{
				int x1 = (1+x0);
				if(x1<nspd)
				{
					bondMap[this->bondMapOffset*bondnum] = x0;
					bondMap[this->bondMapOffset*bondnum + 1] = x1;
					bondMap[this->bondMapOffset*bondnum + 2] = (x0 == nspd-1) ? 1 : 0;

					pairSitesToBond.insert(pair<pair<int,int>,int>(pair<int,int>(bondMap[this->bondMapOffset*bondnum],bondMap[this->bondMapOffset*bondnum+1]),bondnum));
					bondnum++;
				}
				p++;
			}
		}
		else if(this->dimension == 2)
		{
			for(int y0=0;y0<nspd;y0++)
				for(int x0=0;x0<nspd;x0++)
				{
					int x1 = (1 + x0);
					if(x1<nspd)
					{
						bondMap[this->bondMapOffset*bondnum] = p;
						bondMap[this->bondMapOffset*bondnum + 1] = x1 + y0*nspd;
						bondMap[this->bondMapOffset*bondnum + 2] = (x0 == nspd-1) ? 1 : 0;

						pairSitesToBond.insert(pair<pair<int,int>,int>(pair<int,int>(bondMap[this->bondMapOffset*bondnum],bondMap[this->bondMapOffset*bondnum+1]),bondnum));
						bondnum++;
					}

					int y1 = (1 + y0);
					if(y1<nspd)
					{
						bondMap[this->bondMapOffset*bondnum] = p;
						bondMap[this->bondMapOffset*bondnum + 1] = x0 + (y1)*nspd;
						bondMap[this->bondMapOffset*bondnum + 2] = (y0 == nspd-1) ? 2 : 0;

						pairSitesToBond.insert(pair<pair<int,int>,int>(pair<int,int>(bondMap[this->bondMapOffset*bondnum],bondMap[this->bondMapOffset*bondnum+1]),bondnum));
						bondnum++;
					}
					p++; //Go to next site //This can be replaced with x0 + y0*nspd btw
				}
		}
		else if(this->dimension == 3)
		{
			for(int z0=0;z0<nspd;z0++)
				for(int y0=0;y0<nspd;y0++)
					for(int x0=0;x0<nspd;x0++)
					{
						int x1 = (1 + x0);
						if(x1<nspd)
						{

							bondMap[this->bondMapOffset*bondnum] = p;
							bondMap[this->bondMapOffset*bondnum + 1] = x1 + (y0)*nspd + (z0)*nspd*nspd;
							bondMap[this->bondMapOffset*bondnum + 2] = (x1 == nspd-1) ? 1 : 0;

							pairSitesToBond.insert(pair<pair<int,int>,int>(pair<int,int>(bondMap[this->bondMapOffset*bondnum],bondMap[this->bondMapOffset*bondnum+1]),bondnum));
							bondnum++;
						}

						int y1 = (1 + y0);
						if(y1<nspd)
						{
							bondMap[this->bondMapOffset*bondnum] = p;
							bondMap[this->bondMapOffset*bondnum + 1] = x0 + (y1)*nspd + (z0)*nspd*nspd;
							bondMap[this->bondMapOffset*bondnum + 2] = (y1 == nspd-1) ? 2 : 0;

							pairSitesToBond.insert(pair<pair<int,int>,int>(pair<int,int>(bondMap[this->bondMapOffset*bondnum],bondMap[this->bondMapOffset*bondnum+1]),bondnum));
							bondnum++;
						}

						int z1 = (1 + z0);
						if(z1<nspd)
						{
							bondMap[this->bondMapOffset*bondnum] = p;
							bondMap[this->bondMapOffset*bondnum + 1] = x0 + (y0)*nspd + (z1)*nspd*nspd;
							bondMap[this->bondMapOffset*bondnum + 2] = (z1 == nspd-1) ? 3 : 0;

							pairSitesToBond.insert(pair<pair<int,int>,int>(pair<int,int>(bondMap[this->bondMapOffset*bondnum],bondMap[this->bondMapOffset*bondnum+1]),bondnum));
							bondnum++;
						}
						p++;
					}
		}

	}

#ifdef TESTLATTICECREATE
#if DEBUG == 2
	FILE* debug = fopen(DEBUGOUT, "a");
	for(int i=1;i<=ns*dim;i++)
	{
		fprintf(debug,"Bond:%ld => %u %u %u\n",i,bondMap[this->bondMapOffset*i],bondMap[this->bondMapOffset*i + 1],bondMap[this->bondMapOffset*i + 2]);
	}
	fclose(debug);
#endif
#endif
}

void Hamiltonian::initialize()
{
	this->latticeCreate();
	this->createSymmetricSet();
	this->calculateConstant();
}
