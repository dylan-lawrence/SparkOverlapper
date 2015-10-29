/*
 * PairedEndRead.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: qy2
 */

#include "PairedEndRead.h"

PairedEndRead::PairedEndRead() {
	// TODO Auto-generated constructor stub
	this->leftRead = NULL;
	this->rightRead = NULL;
	this->orientation = 0;
	this->mergedSequence = "";
	this->sizeToCount=NULL;
	this->sizeToMatrix=NULL;
	this->pairEndAlignmentList = NULL;

}

PairedEndRead::~PairedEndRead() {
	// TODO Auto-generated destructor stub
	this->clearPairEndAlignmentList();
	this->clearCountMatrix();

}

void PairedEndRead::clearPairEndAlignmentList()
{
	if(this->pairEndAlignmentList!=NULL)
	{
		for(int i=0;i<this->pairEndAlignmentList->size();i++)
			delete this->pairEndAlignmentList->at(i);
		this->pairEndAlignmentList->clear();
		delete this->pairEndAlignmentList;
		this->pairEndAlignmentList = NULL;
	}
}

bool PairedEndRead::addToPairAlignmentList(PairedEndQueryAlignment* pairedAlign)
{
	if(this->pairEndAlignmentList==NULL)
		this->pairEndAlignmentList = new vector<PairedEndQueryAlignment*>();
	this->pairEndAlignmentList->push_back(pairedAlign);
	return true;
}
void PairedEndRead::clearCountMatrix()
{
	if(this->sizeToCount!=NULL)
	{
		this->sizeToCount->clear();
		delete this->sizeToCount;
		this->sizeToCount=NULL;
	}
	if(this->sizeToMatrix!=NULL)
	{
		for (std::map<int, CountMatrix*>::iterator it=this->sizeToMatrix->begin(); it!=this->sizeToMatrix->end(); ++it)
			delete it->second;
		this->sizeToMatrix->clear();
		delete this->sizeToMatrix;
		this->sizeToMatrix=NULL;
	}
}

void PairedEndRead::initializeCountMatrix()
{
	if(this->sizeToCount==NULL)
		this->sizeToCount = new map<int, int>();
	if(this->sizeToMatrix==NULL)
		this->sizeToMatrix = new map<int, CountMatrix*>();
}
bool PairedEndRead::insertMatrixToMap(int readDepthCutoff)
{
	this->initializeCountMatrix();
	if(this->pairEndAlignmentList==NULL || this->pairEndAlignmentList->size()<readDepthCutoff)
		return false;
	else
	{
		for(int i =0; i<this->pairEndAlignmentList->size();i++)
		{
			PairedEndQueryAlignment* align = this->pairEndAlignmentList->at(i);
			int interval = align->right_queryStart-align->left_queryEnd-1;
			if(this->sizeToCount->count(interval)>0)
			{
				(*this->sizeToCount)[interval]=(*this->sizeToCount)[interval]+1;
				CountMatrix* matrix = (*this->sizeToMatrix)[interval];
				matrix->addACTGCount(align);
			}
			else
			{
				int colNum = (interval>0)?interval:-interval;
				if(colNum>0)
				{
				this->sizeToCount->insert(std::pair<int, int>(interval,1));

				CountMatrix * newmatrix = new CountMatrix(5,colNum);
				newmatrix->addACTGCount(align);
				this->sizeToMatrix->insert(std::pair<int, CountMatrix*>(interval, newmatrix));
				}
			}
			delete align;
		}
		this->pairEndAlignmentList->clear();
		delete this->pairEndAlignmentList;
		this->pairEndAlignmentList=NULL;
		return true;
	}


}


bool PairedEndRead::reconstructSequence()
{
	string consensus="";
	int maxCount=0;
	int maxInterval;
	if(this->sizeToCount!=NULL&&this->sizeToCount->size()>0)
	{
		int totalCounts=0;
		for(std::map<int, int>::iterator it=this->sizeToCount->begin(); it!=this->sizeToCount->end(); ++it)
		{
			int curCount = it->second; int index = it->first;
			if(curCount>maxCount)
				{
				totalCounts+=curCount;
				maxCount = curCount;
				maxInterval = index;
				}
		}

		if((double)maxCount/(double)totalCounts<0.8)return false;
		CountMatrix * matrix = (*this->sizeToMatrix)[maxInterval];
		std::stringstream sstm;

		for(int i = 0; i<matrix->matrixColumnSize;i++)
		{
			UINT16 A=matrix->matrixdata->at(0)->at(i);
			UINT16 C=matrix->matrixdata->at(1)->at(i);
			UINT16 T=matrix->matrixdata->at(2)->at(i);
			UINT16 G=matrix->matrixdata->at(3)->at(i);
			UINT16 total = matrix->matrixdata->at(4)->at(i);
			UINT16 maxnum; string maxstring;
			maxnum = A; maxstring = "A";
			if(C>maxnum)
			{
				maxnum = C; maxstring = "C";
			}
			if(T>maxnum)
			{
				maxnum = T; maxstring = "T";
			}
			if(G>maxnum)
			{
				maxnum = G; maxstring = "G";
			}
			if((float)maxnum/(float)total>=0.8&&total>5)sstm<<maxstring;
			else
			{

				return false;
			}
		}
		consensus=sstm.str();
		string left, right="";
		if(maxInterval>0)
		{
			left = this->leftRead->getSequence();
			right = this->rightRead->getSequence();
		}
		else if(maxInterval<0)
		{
			string leftread = this->leftRead->getSequence();
			left = leftread.substr(0,leftread.length()+maxInterval);
			string rightread = this->rightRead->getSequence();
			right = rightread.substr(-maxInterval, rightread.length()+maxInterval);
		}
		this->mergedSequence = left+consensus+right;
	}


	this->clearCountMatrix();
	return true;
}

PairedEndQueryAlignment::PairedEndQueryAlignment(PairedEndRead* pairEndRead,SubjectRead* subjectRead)
{
	this->pairEndRead = pairEndRead;
	this->subjectRead = subjectRead;
	this->left_queryEnd = 0;
	this->right_queryEnd = 0;
	this->right_queryStart = 0;
	this->subjectEnd = 0;
	this->subjectReadOrientation = true;
	this->subjectStart = 0;
}
PairedEndQueryAlignment::~PairedEndQueryAlignment()
{

}

CountMatrix::CountMatrix(int rowSize, int ColumnSize)
{
	this->matrixRowSize = rowSize;
	this->matrixColumnSize = ColumnSize;
	this->matrixdata = new vector<vector<UINT16>*>();
	for(int i= 0; i<this->matrixRowSize; i++)
	{
		vector<UINT16>* matrixRow = new vector<UINT16>();
		for(int j= 0; j<this->matrixColumnSize; j++)
			matrixRow->push_back(0);
		this->matrixdata->push_back(matrixRow);
	}
}

CountMatrix::~CountMatrix()
{
	for(int i= 0; i<this->matrixRowSize; i++)
	{
		vector<UINT16>* matrixRow = this->matrixdata->at(i);
		matrixRow->clear();
		delete matrixRow;
	}
	this->matrixdata->clear();
	delete this->matrixdata;
}

bool CountMatrix::addACTGCount(string sequence)
{
	UINT64 length = sequence.length();// initial value of 1 to avoid zero when doing the multiplication.
	for(UINT64 i = 0; i < length; i++)			// We take the bit representation of the string. A = 00, C = 01, G = 11 and T = 10
	{

			// sum1 is for the first 32 bp. bit shifted to left 2 bits.
			// Change the character to integer. A=Ox41=01000001
			//                                  C=0x43=01000011
			//                                  G=0x47=01000111
			//                                  T=0x54=01010100
			// Then, shift to right way 1 bit.
			// Then, bit and operation with 00000011
			// Then, it just have A=00, C=01, G=11,T=10
		UINT16 rowNumber = (UINT16)(((int)(sequence[i] ) >> 1 ) & 0X03);
		UINT16 colNumber = i;
		vector<UINT16>* currentrow = this->matrixdata->at(rowNumber);
		vector<UINT16>* lastrow = this->matrixdata->at(this->matrixRowSize-1);
		currentrow->at(i) = currentrow->at(i)+1;
		lastrow->at(i) = lastrow->at(i)+1;

	}
	return true;
}
bool CountMatrix::addACTGCount(PairedEndQueryAlignment* pairedEndQueryAlignment)
{
	string subjectSubString;
	PairedEndQueryAlignment* align = pairedEndQueryAlignment;
	int interval = align->right_queryStart-align->left_queryEnd-1;
	int start, end;
	if(interval<0)
	{
		start = align->right_queryStart;
		end = align->left_queryEnd;
	}
	else
	{
		start = align->left_queryEnd+1;
		end = align->right_queryStart-1;
	}

	if(pairedEndQueryAlignment->subjectReadOrientation==true)
	{
		subjectSubString = align->subjectRead->getSequence().substr(start-align->subjectStart, end-start+1);

	}
	else
	{
		string temp = align->subjectRead->getSequence().substr(align->subjectEnd-end, end-start+1);
		subjectSubString = align->subjectRead->reverseComplement(temp);
	}
	return this->addACTGCount(subjectSubString);
}
