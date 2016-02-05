/*
 * DerivedGlobalObservables.h
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

#ifndef DERIVEDGLOBALOBSERVABLES_H_
#define DERIVEDGLOBALOBSERVABLES_H_

namespace measures
{

	class DerivedGlobalObservable : public Measurable
	{
	protected:
		int sample;
		unsigned long opCountavg;
		unsigned long opCount2avg;
		unsigned long densityavg;
		unsigned long density2avg;

	public:

		DerivedGlobalObservable(std::string s, core::StateVariable* lstate, core::Hamiltonian* lham, FILE* llog);
		~DerivedGlobalObservable();

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

#endif /* DERIVEDGLOBALOBSERVABLES_H_ */
