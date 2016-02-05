/*
 * BoseHubbardClean.h
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

#ifndef BOSEHUBBARDCLEAN_H_
#define BOSEHUBBARDCLEAN_H_

class BoseHubbardClean: public core::Hamiltonian {

private:
	double t;
	double U;
	FILE* log;

	double multA;

public:

	BoseHubbardClean(std::string a,FILE*);
	~BoseHubbardClean();

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

#endif /* BOSEHUBBARDCLEAN_H_ */
