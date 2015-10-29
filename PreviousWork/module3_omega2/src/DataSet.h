#ifndef DATASET_H
#define DATASET_H

//============================================================================
// Name        : DataSet.h
// Author      : Tae-Hyuk (Ted) Ahn
// Version     : V0.0.1
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : DataSet header file
//============================================================================

#include "Config.h"
#include "Read.h"

class DataSet
{
	private:
		vector<Read *> *reads;			// List of reads in the dataset.
		void sortReads(void);			// Ted: sort the reads lexicographically.

	public:
		DataSet(void);					// Default constructor.
		~DataSet();						// Destructor.

		UINT64 numberOfUniqueReads; 	// number of unique reads in the dataset.

		// add read into dataset
		void addRead(Read *r){reads->push_back(r);}

		// Get the number of unique reads in the dataset.
		UINT64 getNumberOfUniqueReads(void);

		// Find a read in the database given the ID in constant time.
		Read * getReadFromID(UINT64 ID);
};



#endif /* DATASET_H */
