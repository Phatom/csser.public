/*
 * mainincs.h
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

#ifndef MAININCS_H_
#define MAININCS_H_

	//Basic input/out
	#include <iostream>
	#include <iomanip>
	#include <fstream>
	#include <cstdio>

	//Libraries
	#include <cmath>
	#include <cstdlib>
	#include <ctime>
	#include <cstring>
	#include <gsl/gsl_rng.h>
	#include <gsl/gsl_cdf.h>
	#include <string>
	#include <sstream>
	#include "nr.h"

	//Advanced Data Structures
	#include <unordered_map>
	#include <set>
	#include <list>
	#include <vector>
	#include <map>

	//Sytem libraries
	#include <sys/mman.h>

	//Parallel Libraries
	#include "omp.h"
	#include "mpi.h"

	//Local Libraries
	#include "globalFlags.h"
	#include "constants.h"
	#include "basicDataTypes.h"
	#include "Flags.h"

#endif /* MAININCS_H_ */
