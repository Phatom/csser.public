/*
 * Hamiltonian.h
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

#ifndef HAMILTONIAN_H_
#define HAMILTONIAN_H_

namespace core
{
	struct cmpr
	{
		bool operator()(const std::pair<int, int> a, const std::pair<int, int> b) {
			return (a.first < b.first
					|| (a.first == b.first && a.second < b.second));
		}
	};

	class Hamiltonian
	{
	public:
		unsigned int LAM; //Maximum state flag at site
		double C;
		double EPS;
		double mu;
		int BT; //Types of terms in Hamiltonian

		int dimension;
		long nspd;
		bConditions bc;

		int numSites;
		int numBonds;
		int symmetricSetSize;

		int bondMapOffset;
		unsigned int *bondMap;

		int symmSetOffset;
		unsigned int *symmetricSet; //Contains the sites that represent symmetric groups
		unsigned int *bondToSSet; //Need  something that maps bonds to sites in symmetric group
		std::map<std::pair<int, int>, int, cmpr> pairSitesToBond; //Map site to bond

		std::string latticeParamFileName;


		void latticeCreate();
	public:

		void initialize();

		virtual double operateWith(char op,char state) = 0;
		virtual double getDiagonalTerm(int si,int sj,int statei, int statej) = 0;
		virtual double getOffdiagonalTerm(unsigned int bond, char op, const char* state) = 0;

		virtual int loadHamiltonianParams() = 0;
		virtual void display() = 0;

	protected:
		//The order in which these are called IS VERY IMPORTANT
		virtual int createSymmetricSet() = 0;
		virtual void calculateConstant() = 0;

	};

}
#endif /* HAMILTONIAN_H_ */
