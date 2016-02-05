/*
 * SfluidRotationalObservable.h
 *
 *  Created on: Aug 7, 2014
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

#ifndef SFLUIDROTATIONALOBSERVABLE_H_
#define SFLUIDROTATIONALOBSERVABLE_H_

#include "Measurable.h"

namespace measures {

/*
 *
 */
class SfluidRotationalObservable: public measures::Measurable {
public:

	int size;
	int* siteMap;

	double* rho_s;
	double* Atotal_r;
	double* Atotal;
	int Zcount;

	SfluidRotationalObservable(std::string s, core::StateVariable* lstate, core::Hamiltonian* lham, FILE* llog);
	~SfluidRotationalObservable();

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

} /* namespace measures */
#endif /* SFLUIDROTATIONALOBSERVABLE_H_ */
