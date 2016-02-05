/*
 * stateVariable.h
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

#ifndef STATEVARIABLE_H_
#define STATEVARIABLE_H_

namespace core{

	typedef unsigned int bigi;

	class StateVariable
	{
	public:
		long opCutOff;
		long nOp;

		int statesize;
		char* state;

		long stringsize;
		bigi *opString;
		bigi *vtx;
		bigi *link;

		float beta;
		int lpContrIter;
		bool allocated;

		string stateFileName;
		FILE* stdo;

		float ibeta;

		//Auxilliary variables that are useful
		std::vector<unsigned int> offdiagOperatorLocs;

	public:
		StateVariable();
		StateVariable(int _size,int _strsize, long _opCutOff, FILE* _stdo, float _beta, int _lpc, string _stateFileName);
		~StateVariable();

		virtual int initializeState()=0;

		void reset();
		void allocate();
		int reSize(long newSize);

		int readState(std::string fname);
		int writeState(std::string fname);
		int readState();
		int writeState();

		void display(FILE*);
	};

}
#endif /* STATEVARIABLE_H_ */
