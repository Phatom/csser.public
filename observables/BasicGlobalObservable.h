/*
 * BasicObservable.h
 *
 *  Created on: Apr 22, 2014
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

#ifndef BASICOBSERVABLE_H_
#define BASICOBSERVABLE_H_

namespace measures
{
	struct staticObservable
	{
		unsigned int operatorCount;
		unsigned int density;
		double energy;

		staticObservable() {}

		staticObservable(unsigned int op, unsigned int rho, double e)
		{
			operatorCount = op;
			density = rho;
			energy = e;
		}

	};

	class BasicGlobalObservable : public Measurable
	{
	public:

		double Cbond;

		staticObservable staticO;

		std::vector<staticObservable> collection;

	public:
		 BasicGlobalObservable(std::string s, core::StateVariable* lstate, core::Hamiltonian* lham, FILE* llog);
		 ~BasicGlobalObservable();

		 void measure();
		 void clear();
		 void display();
		 int write();
		 int read();
		 int writeViaIndex(int idx);
		 int readViaIndex(int idx);
		 int checkPointWrite();
		 int checkPointRead();
		 void setup(std::string s, core::StateVariable* lstate, core::Hamiltonian* lham, FILE* llog);
	};
}
#endif /* BASICOBSERVABLE_H_ */
