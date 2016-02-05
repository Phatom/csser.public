/*
 * BHspdmSimEngine.h
 *
 *  Created on: Apr 30, 2014
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

#ifndef BHSPDMSIMENGINE_H_
#define BHSPDMSIMENGINE_H_

namespace measures {

class BHspdmSimEngine: public core::SPDMSimEngine
{

public:
	BHspdmSimEngine(core::Hamiltonian* h, core::TransitionMap* tmap, core::StateVariable* s, gsl_rng* lgsl,
			FILE* llog, map2d& spdmObservable): SPDMSimEngine(h,tmap,s,lgsl,llog,spdmObservable) {}

	int sweep(int&,long long&);
};

} /* namespace measures */
#endif /* BHSPDMSIMENGINE_H_ */
