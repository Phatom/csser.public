/*
 * SSTransitionMap.h
 *
 *	Transition Map for single species of particles
 *	Note that for multiple species the number of ladder operators will change
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

#ifndef SSTRANSITIONMAP_H_
#define SSTRANSITIONMAP_H_

namespace core{

	class SSTransitionMap : public TransitionMap
	{

	public:
		SSTransitionMap(FILE* log, Hamiltonian* ham) : TransitionMap(log,ham) {}
		~SSTransitionMap();

	protected:
		void initializeMaps();
		double inline calcExP(int bond, int si, int sj, long vtxState,
							char ladderOp, char entranceLeg, double* tmatrix);
		void calcExitProbs(std::string base);
	};

}

#endif /* SSTRANSITIONMAP_H_ */
