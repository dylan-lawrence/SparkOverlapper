/*
 * SubjectRead.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "QueryRead.h"

SubjectRead::SubjectRead() {

	readSequence = "";
	readName = "";

	flag4Removal = false;

}

SubjectRead::SubjectRead(string& sequence, string& name)
{
	readSequence = sequence;
	readName = name;

	flag4Removal = false;

}

SubjectRead::~SubjectRead() {
	// TODO Auto-generated destructor stub

}

bool SubjectRead::needRemoval()
{
	return flag4Removal;
}
void SubjectRead::setName(string & name)
{
	readName = name;
}
void SubjectRead::setSequence(string & sequence)
{
	readSequence = sequence;
}

string SubjectRead::getName()
{
	return readName;
}
string SubjectRead::getSequence()
{
	return readSequence;
}


UINT32 SubjectRead::getReadLength()
{
	return this->readSequence.length();
}


string SubjectRead::reverseComplement(const string & read)
{
	UINT64 stringLength = read.length();
	string reverse(stringLength,'0');
	for(UINT64 i = 0;i < stringLength; i++)						// Then complement the string. Change A to T, C to G, G to C and T to A.
	{
		if( read[i] & 0X02 ) // C <==> G
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X04;
		else // A <==> T
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X15;
	}
	return reverse; // return the reverse complement as a string
}

char SubjectRead::reverseBaseAt(const string & read, int position)
{
	if(position<0 || position>=read.length()) return 0;
	char reversebase = 0;
	if( read[position] & 0X02 ) // C <==> G
		reversebase = read[position] ^ 0X04;
	else // A <==> T
		reversebase = read[position] ^ 0X15;
	return reversebase;

}

string SubjectRead::reverseComplement()
{
	return QueryRead::reverseComplement(this->readSequence);
}



