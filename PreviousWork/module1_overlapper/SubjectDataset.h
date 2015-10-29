/*
 * SubjectDataset.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef SUBJECTDATASET_H_
#define SUBJECTDATASET_H_

#include "Config.h"
#include "Alignment.h"
#include "SubjectRead.h"
#include "QueryDataset.h"
//#include "AlignmentPairedEnd.h"

enum FileType {FASTA, FASTQ, UNDEFINED};


class SubjectDataset {

	vector<string> subjectFilenameList;

	UINT64 chunkSize;

	UINT64 currentFileIndex; // the index of the file in sFilenameList, which is currently be processinged loadNextChunk
	FileType currentFileType;

	ifstream currentFileStreamer;
	UINT64 currentChunkSize;

	bool outOfDataFlag;
	bool newFileFlag;

	vector<SubjectRead *>* subjectReadList;
public:
	SubjectDataset();
	~SubjectDataset();


	bool setFilenameList(const vector<string> & sFilenames);
	// load the next chunk and populate currentChunk
	// return false if there is no reads in the file.
	bool loadNextChunk(vector<SubjectRead*> * chunkData);
	bool loadNextChunk(vector<SubjectRead*> * chunkData, QueryDataset * querydataset);
	bool loadNextChunkParallel(vector<SubjectRead*> * chunkData, UINT16 numberOfThreads);
	bool loadNextChunkSingle(vector<SubjectRead*> * chunkData);
//	bool loadNextChunk(vector<SubjectAlignment*> & chunkData);
//	bool loadNextChunk(vector<edge*> & chunkData);
//	bool loadNextChunk(vector<SubjectAlignmentPairedEnd*> & chunkData);
	bool addSubjectRead(SubjectRead* subjectRead);

};

#endif /* SUBJECTDATASET_H_ */
