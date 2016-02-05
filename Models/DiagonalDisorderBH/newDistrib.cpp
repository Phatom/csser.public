/*
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

#define MPI

using namespace std;

extern int equibRunner(int rank, string path, char* argv);
extern int warmRunner(int rank, string path, char* argv);
extern int dataRunner(int rank, string path, char* argv);

int main(int argc, char* argv[])
{
	if(argc<3)
	{
		cout<<"Format is: <Mode> <Parameter File> [<Directory List>]"<<endl;
		cout<<"Valid modes are: Eqb, Warmup, Data\n";
		return -1;
	}


#ifdef MPI
	if(argc<4)
	{
		cout<<"Format is: <Mode> <Parameter File> <Directory List>"<<endl;
		return -1;
	}

	int rank;
	int prov;

	MPI_Status Stat;
	MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE, &prov);

	if(prov == MPI_THREAD_SINGLE)
		cout<<"Thread Support: SINGLE"<<endl;
	else if(prov == MPI_THREAD_FUNNELED)
		cout<<"Thread Support: FUNNEL"<<endl;
	else if(prov == MPI_THREAD_SERIALIZED)
		cout<<"Thread Support: SERIALIZED"<<endl;
	else if(prov == MPI_THREAD_MULTIPLE)
		cout<<"Thread Support: MULTIPLE"<<endl;
	else
		cout<<"Thread Support: ERROR"<<endl;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	string prefix = "";
	ifstream pf(argv[3]);
	bool done = false;
	for(int i=0;i<=rank;i++)
	{
		pf>>prefix;
		if(!pf)
		{
			done = true;
			break;
		}
	}
	pf.close();

	if(done)
	{	
		MPI_Finalize();
		return 0;
	}	

	//Setup file first
	cout<<"Prefix is: "<<prefix<<endl;

	//Now call modules as needed
	cout<<"Staring with mode: "<<argv[1]<<endl;
	if(string(argv[1]).compare("Eqb") == 0)
		equibRunner(rank, prefix, argv[2]);
	else if(string(argv[1]).compare("Warmup") == 0)
		warmRunner(rank, prefix, argv[2]);
	else if(string(argv[1]).compare("Data") == 0)
		dataRunner(rank, prefix, argv[2]);

	MPI_Finalize();
#else

	//No MPI
	int rank = 0;
	string prefix = "";

	//Now call modules as needed
	cout<<"Staring with mode: "<<argv[1]<<endl;
	if(string(argv[1]).compare("Eqb") == 0)
		equibRunner(rank, prefix, argv[2]);
	else if(string(argv[1]).compare("Warmup") == 0)
		warmRunner(rank, prefix, argv[2]);
	else if(string(argv[1]).compare("Data") == 0)
		dataRunner(rank, prefix, argv[2]);

#endif

}
