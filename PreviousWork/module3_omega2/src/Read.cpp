//============================================================================
   // Name        : Read.cpp
// Author      : Tae-Hyuk (Ted) Ahn
// Version     : V0.0.1
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Read cpp file
//============================================================================


#include "Config.h"
#include "Read.h"


//=============================================================================
// Default constructor
//=============================================================================
Read::Read(void)
{
	// Initialize the variables.
	readNumber = 0;
	frequency = 0;
	listOfEdgesForward = new vector<Edge *>;
	listOfEdgesForward->resize(listOfEdgesForward->size());			// Resize to 0 to reduce space.
	listOfEdgesReverse = new vector<Edge *>;
	listOfEdgesReverse->resize(listOfEdgesReverse->size());			// Resize to 0 to reduce space.
	locationOnEdgeForward = new vector<UINT64>;
	locationOnEdgeForward->resize(locationOnEdgeForward->size());	// Resize to 0 to reduce space.
	locationOnEdgeReverse = new vector<UINT64>;
	locationOnEdgeReverse->resize(locationOnEdgeReverse->size());	// Resize to 0 to reduce space.
}


//=============================================================================
// Another constructor
//=============================================================================
Read::Read(const string & s)
{
	// Initialize the variables.
	readNumber = 0;
	frequency = 0;
	listOfEdgesForward = new vector<Edge *>;
	listOfEdgesForward->resize(listOfEdgesForward->size());			// Resize to 0 to reduce space.
	listOfEdgesReverse = new vector<Edge *>;
	listOfEdgesReverse->resize(listOfEdgesReverse->size());			// Resize to 0 to reduce space.
	locationOnEdgeForward = new vector<UINT64>;
	locationOnEdgeForward->resize(locationOnEdgeForward->size());	// Resize to 0 to reduce space.
	locationOnEdgeReverse = new vector<UINT64>;
	locationOnEdgeReverse->resize(locationOnEdgeReverse->size());	// Resize to 0 to reduce space.
	setRead(s);
}


//=============================================================================
// Default destructor
//=============================================================================
Read::~Read(void)
{
	// delete all the pointers.
	delete listOfEdgesForward;
	delete listOfEdgesReverse;
	delete locationOnEdgeForward;
	delete locationOnEdgeReverse;
}


//=============================================================================
//	Function to store the read sequence
//=============================================================================
bool Read::setRead(const string & s)
{
	setFrequency(1);												// Set the frequency to 1.
	read = s;														// Store the string.
	readReverse = reverseComplement(s);
	return true;
}


//=============================================================================
// Function to store the read number (ID)
//=============================================================================
bool Read::setReadNumber(UINT64 id)
{
	if(id <= 0) MYEXIT("ID less than 1.");
	readNumber = id;												// Set the read number.
	return true;
}


//=============================================================================
// Function to set the frequency of the read.
//=============================================================================
bool Read::setFrequency(UINT32 freq)
{
	if(freq < 1) MYEXIT("Frequency less than 1.");
	frequency = freq;												// Set the frequency of the read.
	return true;
}


//=============================================================================
// Returns the reverse complement of a read.
//=============================================================================
string Read::reverseComplement(const string & read)
{
	UINT64 stringLength = read.length();
	string reverse(stringLength,'0');
	for(UINT64 i = 0;i < stringLength; i++)						// Then complement the string. Change A to T, C to G, G to C and T to A.
	{
		if( read[i] & 0X02 ) // C or G
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X04;
		else // A <==> T
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X15;
	}
	return reverse; // return the reverse complement as a string
}
