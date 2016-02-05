/*
 * Observable.h
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

#ifndef MEASURABLE_H_
#define MEASURABLE_H_

namespace measures
{
	class Measurable
	{
	protected:
		std::string fileName, chkFileName;
		core::StateVariable *state;
		core::Hamiltonian *hamiltonian;

		FILE* log;

	public:
		Measurable(std::string s, core::StateVariable* lstate, core::Hamiltonian* lham, FILE* llog);

		void setState(core::StateVariable *s);
		void setFileName(std::string fname);
		void setLog(FILE* llog);

		virtual void measure() = 0;
		virtual void clear() = 0;
		virtual void display() = 0;

		virtual int writeViaIndex(int idx) = 0;
		virtual int readViaIndex(int idx) = 0;
		virtual int write() = 0;
		virtual int read() = 0;

		virtual int checkPointWrite() = 0;
		virtual int checkPointRead() = 0;
	};
}
#endif /* OBSERVABLE_H_ */
