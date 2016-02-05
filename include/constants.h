/*
 * constants.h
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

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

///////////////////////////////////////
// Size Cutoffs
///////////////////////////////////////

#define HARDCUTOFF (1.5e9)

//Basic error codes
#define EXITCODE -3
#define FILENOTFOUND -2
#define PARTIALFAIL -1
#define SUCCESS 0
#define LOOPLENGTHEXCEEDED -4
#define ConstructionException -5
//Float Precsions and Widths
#define FIELDPRECISION 9
#define FIELDWIDTH 10
#define FIELDFORMAT std::ios_base::scientific
#endif /* CONSTANTS_H_ */
