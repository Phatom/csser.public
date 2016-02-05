/*
 * flagEnums.h
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

#ifndef FLAGENUMS_H_
#define FLAGENUMS_H_

enum bConditions {HW, PBC};
///////////////////////////////////////////////////////////////////////

struct eqnum{
	bool operator()(int num1, int num2)
	const
	{
		return num1 == num2;
	}
};

struct gcmpr{
	bool operator()(const int a, const int b)
	{
		return a<b;
	}
};

typedef unordered_map<int,double,hash<int>,eqnum> Row;
typedef map<long,double,gcmpr> sortedRow;

typedef unordered_map<int,sortedRow,hash<int>,eqnum> orderedMap2d;
typedef unordered_map<int,Row,hash<int>,eqnum> hashmap2d;

typedef unordered_map<int,hashmap2d,hash<int>,eqnum> hashmap3d;

typedef hashmap2d map2d;
typedef hashmap3d map3d;


struct ival
{
	int endLevel;
	float value;

	ival(): endLevel(0), value(0.0) {}
	ival(int l, float v): endLevel(l), value(v) { }
	ival(const ival& obj)
	{
		endLevel = obj.endLevel;
		value = obj.value;
	}
};

typedef unordered_map<int,ival,hash<int>,eqnum> Interval;
typedef unordered_map<int,Interval,hash<int>, eqnum>  SiteToInterval;
typedef unordered_map<int,SiteToInterval,hash<int>,eqnum> PairSiteToInterval;
typedef PairSiteToInterval GijI;

#endif /* FLAGENUMS_H_ */
