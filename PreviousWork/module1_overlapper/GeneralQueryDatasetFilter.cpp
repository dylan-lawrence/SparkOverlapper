/*
 * GeneralQueryDatasetFilter.cpp
 *
 *  Created on: Apr 29, 2015
 *      Author: qy2
 */

#include "GeneralQueryDatasetFilter.h"

GeneralQueryDatasetFilter::GeneralQueryDatasetFilter(HashTableMethod * hashTableMethod) {
	// TODO Auto-generated constructor stub

	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashTableMethod = hashTableMethod;
}

GeneralQueryDatasetFilter::~GeneralQueryDatasetFilter() {
	// TODO Auto-generated destructor stub
}

bool GeneralQueryDatasetFilter::start()
{
	string outFileName = Config::outputfilename;
	filePointer.open(outFileName.c_str());
	CLOCKSTART;
	MEMORYSTART;

	vector<SubjectRead *>* subjectReadList= new vector<SubjectRead *>();
	subjectReadList->clear();

	SubjectDataset *subjectDataset = new SubjectDataset();
	subjectDataset->setFilenameList(Config::subjectFilenameList);
	while(subjectDataset->loadNextChunk(subjectReadList))
	{
		vector<SubjectEdge*>* subjectEdgeList = new  vector<SubjectEdge*>();
		for(UINT64 i = 0; i < subjectReadList->size(); i++)
		{
			SubjectRead * subjectRead = subjectReadList->at(i);
			SubjectEdge * subjectEdge = new SubjectEdge(subjectRead);
			subjectEdgeList->push_back(subjectEdge);
		}
omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(UINT64 i = 0; i < subjectEdgeList->size(); i++)
		{
			SubjectEdge * subjectEdge = subjectEdgeList->at(i);
			this->hashTableMethod->searchHashTableFullScanAndMarkQueryRead(subjectEdge); //and also mark the contained read in the query read structure
			delete subjectEdge->subjectRead;
			delete subjectEdge;
		}

	}
		// end of omp parallel

		//
		subjectEdgeList->clear();
		delete subjectEdgeList;
		//clear the subject read vector for the next round
		subjectReadList->clear();
		subjectReadList->resize(0);

	}// end of while chunk loop
	this->printToFile();

	delete subjectReadList;
	delete subjectDataset;

	MEMORYSTOP;
	CLOCKSTOP;

filePointer.close();
	return true;
}

void GeneralQueryDatasetFilter::printToFile()
{
	UINT64 count = 0;
	QueryDataset * dataset = this->hashTableMethod->getDataset();
	if(Config::useID == false)
	{
		for(UINT64 i = 1; i <= dataset->getNumberOfUniqueReads(); i++)
		{
			QueryRead *queryRead = dataset->getReadFromID(i);
			if(queryRead->flag4Removal==false)
			{
				count++;
				std::stringstream sstm;
				sstm<<">"<<queryRead->getName()<<endl;
				sstm<<queryRead->getSequence()<<endl;
				this->filePointer<<sstm.str();
				if(count%100000==0)	// Show the progress.
					cout<<"counter: " << setw(10) << count << " Reads Written. " << setw(10) << endl;
			}
		}
	}
	else
	{
		for(UINT64 i = 1; i <= dataset->getNumberOfUniqueReads(); i++)
		{
			QueryRead *queryRead = dataset->getReadFromID(i);
			if(queryRead->flag4Removal==false)
			{
				count++;
				std::stringstream sstm;
				UINT64 queryReadID = queryRead->getIdentifier()+Config::startID;
				sstm<<">"<<queryReadID<<endl;
				sstm<<queryRead->getSequence()<<endl;
				this->filePointer<<sstm.str();
				if(count%100000==0)	// Show the progress.
					cout<<"counter: " << setw(10) << count << " Reads Written. " << setw(10) << endl;
			}
		}
	}
	cout<<"Total counter: " << setw(10) << count << " Reads Written. " << setw(10) << endl;
}
