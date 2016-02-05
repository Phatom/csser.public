/*
 * BasicLocalObservable.h
 *
 *  Created on: Apr 23, 2014
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

#ifndef BASICLOCALOBSERVABLE_H_
#define BASICLOCALOBSERVABLE_H_

namespace measures
{
	class BasicLocalObservable : public Measurable
	{
	protected:
		int bin;
		int sample;
		float* localRho;
		float* localRho2;

	public:
		BasicLocalObservable(std::string s, core::StateVariable* lstate, core::Hamiltonian* lham, FILE* llog);
		~BasicLocalObservable();

		 void measure();
		 void clear();
		 void display();
		 int write();
		 int read();
		 int writeViaIndex(int idx);
		 int readViaIndex(int idx);
		 int checkPointWrite();
		 int checkPointRead();
	};
}


#endif /* BASICLOCALOBSERVABLE_H_ */
