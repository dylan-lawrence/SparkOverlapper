//============================================================================
// Name        : DataSet.cpp
// Author      : Tae-Hyuk (Ted) Ahn
// Version     : V0.0.1
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : DataSet cpp file
//============================================================================

#include "Config.h"
#include "DataSet.h"



//=============================================================================
// Default constructor
//=============================================================================
DataSet::DataSet(void)
{
	reads = new vector<Read *>;
	numberOfUniqueReads = 0;
}

//=============================================================================
// Default destructor
//=============================================================================
DataSet::~DataSet()
{
	// Free the memory used by the dataset.
	for(UINT64 i = 0; i < reads->size(); i++)
	{
		delete reads->at(i);
	}
	delete reads;
}

//=============================================================================
// This function returns the number of unique reads
//=============================================================================
UINT64 DataSet::getNumberOfUniqueReads(void)
{
	return numberOfUniqueReads;
}


//=============================================================================
//	This function returns a read from its ID.
//=============================================================================
Read * DataSet::getReadFromID(UINT64 ID)
{
	if(ID < 1 || ID > numberOfUniqueReads)	// ID outside the range.
	{
		stringstream ss;
		ss << "ID " << ID << " out of bound.";
		string s = ss.str();
		MYEXIT(s);
	}
	return reads->at(ID - 1);
}
