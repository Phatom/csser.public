/*
 * BoseHubbardDDisorder.h
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

#ifndef BoseHubbardDDisorder_H_
#define BoseHubbardDDisorder_H_

class BoseHubbardDDisorder: public core::Hamiltonian {

private:
	int processId;
	string path;

	double t;
	double U;
	double* E;
	double Delta;
	double Omega;
	int disorderSeed;

	FILE* log;
	double multA;
	gsl_rng* rgenref;
	string disorderFileName;

public:

	BoseHubbardDDisorder(int _procId,std::string path,std::string hamFilename,FILE* log);
	~BoseHubbardDDisorder();

	void setLatticeFileName(string fn);

	double operateWith(char op,char state) ;
	double getDiagonalTerm(int si,int sj,int statei, int statej) ;
	double getOffdiagonalTerm(unsigned int bond, char op, const char* state) ;

	int loadHamiltonianParams();
	void display();

protected:
	int createSymmetricSet();
	void calculateConstant();

};

#endif /* BoseHubbardDDisorder_H_ */
