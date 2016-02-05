/*
 * GreenObservable.h
 *
 *  Created on: May 3, 2014
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



#ifndef GREENDYNAMEASURABLE_H_
#define GREENDYNAMEASURABLE_H_

namespace measures {

class GreenObservable: public measures::Measurable {

protected:

	Flags flag;

	orderedMap2d tvalues;
	int numTvalues;
	vector<float> tcollection;

	float sigWidth;
	long L;

	map3d green;
	GijI &lpvar;

	int Zcount;



public:
	 GreenObservable(std::string s, core::StateVariable* lstate, core::Hamiltonian* lham, FILE* llog, Flags& f, GijI& _loopVar);
	 ~GreenObservable();

	 void createTvalues();
	 void measure();
	 void clear();
	 void display();
	 int write();
	 int read();
	 int writeViaIndex(int);
	 int readViaIndex(int);
	 int checkPointWrite();
	 int checkPointRead();
};

} /* namespace measures */
#endif /* SPDMDYNAMEASURABLE_H_ */
