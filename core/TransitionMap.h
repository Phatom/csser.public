/*
 * TransitionMap.h
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

#ifndef TRANSITIONMAP_H_
#define TRANSITIONMAP_H_

namespace core {

	class TransitionMap {
		protected:

			bool mapsInitialized;
			char* opmap;//This is about what is the operator given 1) entrance leg 2) exit leg 3) current operator type
			char** statemap;	//Given vertex type what are the states of the four legs
			long* vtxTypeMap;	//Given  state of 4 legs what is vertex type
			long* newVtxMap;//Given vertex type, entrance leg no. and exit leg no., and operator what is new vertex

			bool probsCalculated;
			double* bondToProb;	//Given bond what are the vertex transitions


			//pointer shifts
			int tsize;			//Total size of maps
			int vtxCount; 		//Keeps track of the number of valid vertices

			int VM, VM2, VM3, VM4; //For tracking state offsets
			int vtxOffset;
			int symmOffset;	//Track symmetries for the general case = number of bonds i.e. all bonds are different

			Hamiltonian* hamiltonian; //Reference to Hamiltonian
			FILE* log; //output stream

			///////////////////////////////////////////////////////////////////////////////
			friend class Hamiltonian;
			friend class SimulationEngine;
			friend class BasicSimEngine;
			friend class SPDMSimEngine;
			friend class GreenSimEngine;
			///////////////////////////////////////////////////////////////////////////////

		public:

			TransitionMap(FILE* _log, Hamiltonian* _ham)
			{
				log = _log;
				hamiltonian = _ham;
				tsize = 0;
				VM=VM2=VM3=VM4=0;
				vtxCount = 0;
				vtxOffset = symmOffset = 0;

				probsCalculated = false;
				mapsInitialized = false;
			}

			~TransitionMap() {}

			void writeBondProbs(std::string base);
			int readBondProbs(std::string base);
			void initialize(std::string base);

		protected:
			virtual void initializeMaps() = 0;
			virtual double inline calcExP(int bond, int si, int sj, long vtxState,
					char ladderOp, char entranceLeg, double* tmatrix) = 0;
			virtual void calcExitProbs(std::string base) = 0;
	};
}

#endif /* TRANSITIONMAP_H_ */
