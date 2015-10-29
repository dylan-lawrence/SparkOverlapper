#ifndef READ_H
#define READ_H

//============================================================================
// Name        : Read.h
// Author      : Tae-Hyuk (Ted) Ahn
// Version     : V0.0.1
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Read header file
//============================================================================

#include "Config.h"

class Edge;


class Read
{
	private:
		UINT64 readNumber; 						// Unique Identification of the read.
		string read; 							// String representation of the read.
		string readReverse;
		UINT32 frequency; 						// Frequency of the read. Number of times this read is present in the dataset. Used in some statistical analysis.
		// listOfEdgesForward and locationOnEdgeReverse stores the list of forward version and reverse version of the same set of edges
		vector<Edge *> *listOfEdgesForward;   	// List of edges that contain the forward string of this read.
		// Location is the distance from the source read on the edge to the beginning of this read.
		vector<UINT64> *locationOnEdgeForward;	// List of locations on the edges that contain the forward string of the current read.
		vector<Edge *> *listOfEdgesReverse;		// List of edges that contain the reverse string of this read.
		vector<UINT64> *locationOnEdgeReverse;	// List of locations on the edges that contain the reverse string of the current read.

		string reverseComplement(const string & read);

	public:
		Read(void);								// Default constructor.
		Read(const string & s);					// Another constructor.
		~Read(void);							// Destructor.

		bool setReadNumber(UINT64 id); 			// Set the read number.
		bool setRead(const string & s); // Set the read sequence.
		bool setFrequency(UINT32 freq);			// Set the frequency of the read.

		string getStringForward(void) const {return read;} 									// Get the forward string of the current read.
		string getStringReverse(void) const {return readReverse;} 								// Get the reverse string of the current read.
		UINT32 getReadLength(void) const {return read.length();} 								// Get the length of the string in the current read.
		UINT64 getReadNumber(void) const {return readNumber;} 								// Get the read number of the current read.
		UINT32 getFrequency(void) const {return frequency;}									// Get the frequency of the current read.
		vector<Edge *> * getListOfEdgesForward(void) const {return listOfEdgesForward;}		// Get the list of edges that contain the forward string of the current read.
		vector<UINT64> * getLocationOnEdgeForward(void) const {return locationOnEdgeForward;}	// Get the list of locations on the edges that contain the forward string of the current read.
		vector<Edge *> * getListOfEdgesReverse(void) const {return listOfEdgesReverse;}		// Get the list of edges that contain the reverse string of the current read.
		vector<UINT64> * getLocationOnEdgeReverse(void) const {return locationOnEdgeReverse;}	// Get the list of locations on the edges that contain the reverse string of the current read.
		bool addMatePair(UINT64 read, UINT8 orientation, UINT64 datasetNumber);				// Add a matepair in the list.

};

#endif /* READS_H */
