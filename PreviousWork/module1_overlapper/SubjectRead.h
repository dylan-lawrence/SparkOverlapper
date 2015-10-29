/*
 * QueryRead.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef SUBJECTREAD_H_
#define SUBJECTREAD_H_

#include "Config.h"
#include "Alignment.h"
class Alignment;
class SubjectRead {
	string readSequence;
	string readName;


public:
	bool flag4Removal; // this is a contained or duplicate read
	SubjectRead();
	SubjectRead(string & sequence, string & name);
	~SubjectRead();
	bool needRemoval();
	void setName(string & name);
	void setSequence(string & sequence);
	string getName();
	string getSequence();
	UINT32 getReadLength();
	static string reverseComplement(const string & read);
	string reverseComplement();
	char reverseBaseAt(const string & read, int position);
};

#endif /* SUBJECTREAD_H_ */
