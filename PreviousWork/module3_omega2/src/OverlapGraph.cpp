//============================================================================
// Name        : OverlapGraph.cpp
// Author      : Tae-Hyuk (Ted) Ahn
// Version     : V0.0.1
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : OverlapGraph cpp file
//============================================================================

#include "Config.h"
#include "OverlapGraph.h"
//#include "CS2/cs2.h"
#include "CS2_stream/cs2.h"

//=============================================================================
// Check if two edges match.
// e1(u,v) and e2(v,w). At node v, one of the edges should be an incoming edge and the other should be an outgoing
// edge to match.
//=============================================================================
bool matchEdgeType(const Edge *edge1, const Edge *edge2)
{
	if     ( (edge1->getOrientation() == 1 || edge1->getOrientation() == 3) && (edge2->getOrientation() == 2 || edge2->getOrientation() == 3) ) // *-----> and >------*
		return true;
	else if( (edge1->getOrientation() == 0 || edge1->getOrientation() == 2) && (edge2->getOrientation() == 0 || edge2->getOrientation() == 1) ) // *------< and <------*
		return true;
	return false;
}


//=============================================================================
// Function to compare two edges. Used for sorting.
//=============================================================================
bool compareEdgeID (const Edge *edge1, const Edge* edge2)
{
	return (edge1->getDestinationRead()->getReadNumber() < edge2->getDestinationRead()->getReadNumber());
}


//=============================================================================
// Function to compare two edges. Used for sorting.
//=============================================================================
// CP: Sort by overlap offset. What for?
// BH: In the transitive edge removal step, we need the edges to be sorted according to their overlap offset.
bool compareEdges (const Edge *edge1, const Edge* edge2)
{
	return (edge1->getOverlapOffset() < edge2->getOverlapOffset());
}

//=============================================================================
// Function to compare two string length
//=============================================================================
bool byStringLength (const string seq1, const string seq2)
{
	return (seq1.length() < seq2.length());
}


//=============================================================================
// Default constructor
//=============================================================================
OverlapGraph::OverlapGraph(void)
{
	// Initialize the variables.
	numberOfNodes = 0;
	numberOfEdges = 0;
	flowComputed = false;
}

//=============================================================================
// Default destructor
//=============================================================================
OverlapGraph::~OverlapGraph()
{
	// Free the memory used by the overlap graph.
	for(UINT64 i = 0; i < graph->size(); i++)
	{
		for(UINT64 j = 0; j< graph->at(i)->size(); j++)
		{
			delete graph->at(i)->at(j);
		}
		delete graph->at(i);
	}
	delete graph;
}


//=============================================================================
// Build overlap graph from read/edge file(s)
//=============================================================================
bool OverlapGraph::buildOverlapGraphFromFiles(const vector<string> &readFilenameList, const vector<string> &edgeFilenameList)
{
	CLOCKSTART;

	// Reads
	UINT64 readID = 0;
	UINT64 numOfUniqueRead = 0;			// Initialize
	UINT64 numOfUniqueReadEachFile = 0;	// Initialize

	// save read ID from file
	//vector<UINT64> readIDListInFile;
	unordered_map<UINT64, UINT64> readIDMap;

	// loop readFilenameList
	for (vector<string>::const_iterator it = readFilenameList.begin(); it != readFilenameList.end(); ++it) {
		readID = numOfUniqueRead+1;
		//numOfUniqueReadEachFile = loadReadFile(*it, readIDListInFile, readID);
		numOfUniqueReadEachFile = loadReadFile(*it, readIDMap, readID);
		numOfUniqueRead += numOfUniqueReadEachFile;
	}

	if (numOfUniqueRead == readIDMap.size()) {
		cout << "numOfUniqueRead=" << numOfUniqueRead << endl;
		cout << "readIDListInFile.size=" << readIDMap.size() << endl;
		dataSet->numberOfUniqueReads = numOfUniqueRead;
	}
	else {
		MYEXIT("Reads load error!");
	}

	// Initialize graph
	graph = new vector< vector<Edge *> * >;
	graph->reserve(dataSet->getNumberOfUniqueReads()+1);
	for(UINT64 i = 0; i <= dataSet->getNumberOfUniqueReads(); i++) // initialize the list
	{
		vector<Edge *> *newList = new vector<Edge *>;
		graph->push_back(newList);
	}

	//removedEdges = new vector <Edge>;
	contigSequences = new vector <string>;

	// loop edgeFilenameList
	for (vector<string>::const_iterator it = edgeFilenameList.begin(); it != edgeFilenameList.end(); ++it) {
    	loadEdgeFile(*it, readIDMap);
	}
	cout << "numberOfEdges = " << numberOfEdges << endl;

	// delete vector
	//readIDListInFile.clear();
	readIDMap.clear();

    // Sort all its incident edges according to destination read ID.
    sortEdges();

    // Composite edge contraction with remove dead end nodes
    UINT64 counter = 0;
	do
	{
		 counter = contractCompositePaths();
		 counter += removeDeadEndNodes();
		 cout << "numberOfEdges = " << numberOfEdges << endl;
	} while (counter > 0);


	CLOCKSTOP;
	return true;
}


//=============================================================================
// Load read file
//=============================================================================
//UINT64 OverlapGraph::loadReadFile(const string &readFilename, vector<UINT64> &readIDListInFile, UINT64 &readID)
UINT64 OverlapGraph::loadReadFile(const string &readFilename, unordered_map<UINT64, UINT64> &readIDMap, UINT64 &readID)
{
	CLOCKSTART;
#ifdef DEBUG
	cout << "readFilename: " << readFilename << endl;
#endif

	// To count of reads
	UINT64 readCount = 0;

	// Open file
    ifstream filePointer;
    filePointer.open(readFilename.c_str());
    if(filePointer == NULL)
        MYEXIT("Unable to open file: "+readFilename);

    // Variables
    vector<string> line;
    string line0,line1, text;
    enum FileType {FASTA, FASTQ, UNDEFINED};
    FileType fileType = UNDEFINED;

    while(!filePointer.eof())
    {
    	// Check FASTA and FASTQ
        if(fileType == UNDEFINED)
        {
            getline (filePointer,text);
            if(text[0] == '>')
                fileType = FASTA;
            else if(text[0] == '@')
                fileType = FASTQ;
            else
                cout<< "Unknown input file format."<<endl;
            filePointer.seekg(0, ios::beg);
        }

        // Clear line
        line.clear();

        // FASTA file read
        if(fileType == FASTA)
        {
        	getline (filePointer,text);		// get ID line
            line.push_back(text);
            getline (filePointer,text,'>');	// get string line
            line.push_back(text);

            line.at(1).erase(std::remove(line.at(1).begin(), line.at(1).end(), '\n'), line.at(1).end());
            line0 = line.at(0);		// ID
            line1 = line.at(1);		// String

        }
        // FASTQ file read
        else if(fileType == FASTQ)
        {
            for(UINT64 i = 0; i < 4; i++)   // Total of 4 lines represent one sequence in a FASTQ file.
            {
                getline (filePointer,text);
                line.push_back(text);
            }
            line0 = line.at(0);		// ID
            line1 = line.at(1);		// String
        }

        // Get ReadID after removing the > or @ identifier and convert string to UINT64
        string readName="";
        if(line0[0] == '>' || line0[0] == '@')
            readName = line0.substr(1);
        else readName = line0;

        istringstream readNameStream(readName);
        UINT64 readIDInFile;
        readNameStream >> readIDInFile;
        //readIDListInFile.push_back(readIDInFile);
        readIDMap[readIDInFile] = readID;

        // Capitalize the sequence
        string read;
        for (std::string::iterator p = line1.begin(); line1.end() != p; ++p)
            *p = toupper(*p);
        read = line1;

        // Save to Read object
        Read *r=new Read();
        // readID, read(sequence), frequency
        r->setReadNumber(readID);
        r->setRead(read);
        r->setFrequency(1);		// Set frequency to 1 for all reads
        // add read to the dataset
        dataSet->addRead(r);
        readCount ++;
        readID ++;
    }
    filePointer.close();
    CLOCKSTOP;

    return readCount;
}


//=============================================================================
// Load edge file
//=============================================================================
void OverlapGraph::loadEdgeFile(const string &edgeFilename, unordered_map<UINT64, UINT64> &readIDMap)
{
	CLOCKSTART;
#ifdef DEBUG
	cout << "edgeFilename: " << edgeFilename << endl;
#endif

	// Open file
    ifstream filePointer;
    filePointer.open(edgeFilename.c_str());
    if(filePointer == NULL)
        MYEXIT("Unable to open file: "+edgeFilename);

    // Read file
    string text;
    while(getline(filePointer,text))
	{
		// Save text to line vector
		vector<string> line;
		std::stringstream text_ss(text);
		std::string element;

		// loop
		while (getline(text_ss, element, '\t'))
		{
			line.push_back(element);
		}

		// Get source, destination, and properties
		UINT64 source_org, source;			// Tab delimited Column1: source read ID
		istringstream t1(line.at(0)); t1 >> source_org;
		//source = std::find(readIDListInFile.begin(), readIDListInFile.end(), source_org) - readIDListInFile.begin() + 1;	// search index and add 1 as readID starts with 1
		unordered_map<UINT64, UINT64>::const_iterator got_s = readIDMap.find (source_org);
		if ( got_s == readIDMap.end() ) {
			MYEXIT("not found");
		}
		else {
			source = got_s->second;
		}

		UINT64 destination_org, destination;		// Tab delimited Column2: destination read ID
		istringstream t2(line.at(1)); t2 >> destination_org;
		//destination = std::find(readIDListInFile.begin(), readIDListInFile.end(), destination_org) - readIDListInFile.begin() + 1;
		unordered_map<UINT64, UINT64>::const_iterator got_d = readIDMap.find (destination_org);
		if ( got_d == readIDMap.end() ) {
			MYEXIT("not found");
		}
		else {
			destination = got_d->second;
		}

		string properties = line.at(2);		// Tab delimited Column3: edge attributes

		// For properties list
		vector<string> propertiesList;
		std::stringstream properties_ss(properties);
		while (getline(properties_ss, element, ','))
		{
			propertiesList.push_back(element);
		}

		// properties
		UINT32 orientation;		// Property Col1: orientation
		std::istringstream p1(propertiesList.at(0)); p1 >> orientation;
		UINT32 overlapLength;	// Property Col2: overlap length
		std::istringstream p2(propertiesList.at(1)); p2 >> overlapLength;
		UINT32 substitutions;	// Property Col3: substitutions
		std::istringstream p3(propertiesList.at(2)); p3 >> substitutions;
		UINT32 edits;			// Property Col4: edit distance
		std::istringstream p4(propertiesList.at(3)); p4 >> edits;
		UINT32 length1;			// Property Col5: read1 length
		std::istringstream p5(propertiesList.at(4)); p5 >> length1;
		INT32 start1;			// Property Col6: read1 overlap start
		std::istringstream p6(propertiesList.at(5)); p6 >> start1;
		INT32 stop1;			// Property Col7: read1 overlap end
		std::istringstream p7(propertiesList.at(6)); p7 >> stop1;
		UINT32 length2;			// Property Col8: read2 length
		std::istringstream p8(propertiesList.at(7)); p8 >> length2;
		INT32 start2;			// Property Col9: read2 overlap start
		std::istringstream p9(propertiesList.at(8)); p9 >> start2;
		INT32 stop2;			// Property Col10: read2 overlap length
		std::istringstream p10(propertiesList.at(9)); p10 >> stop2;
		string errorInfo;		// Property Col11: error info
		errorInfo = propertiesList.at(10);

		// get overlap offset
		// Example edge list:
		// 0 = u<-----------<v		reverse of u to reverse of v
		// 1 = u<----------->v		reverse of u to forward of v
		// 2 = u>-----------<v		forward of u to reverse of v
		// 3 = u>----------->v		forward of u to forward of v
		//1	2   2,34,0,0,35,1,34,35,0,33,NA
		//2	3   0,33,0,0,35,2,34,35,0,32,NA
		//3	4	3,34,1,1,35,1,34,35,0,33,1CA

		// overlap offset
		UINT32 overlapOffset;
		overlapOffset = start1;	// correct, but not ready yet
		// Insert an edge
		Read *read1, *read2;
		read1 = dataSet->getReadFromID(source);
		read2 = dataSet->getReadFromID(destination);
		insertEdge(read1,read2,orientation,overlapOffset);			// insert edge
    }
    filePointer.close();

    CLOCKSTOP;
}


//=============================================================================
// Insert an edge in the overlap graph.
//=============================================================================
bool OverlapGraph::insertEdge(Edge * edge)
{

	UINT64 ID = edge->getSourceRead()->getReadNumber(); // This is the source read.
	if(graph->at(ID)->empty()) 							// If there is no edge incident to the node
		numberOfNodes++;								// Then a new node is inserted in the graph. Number of nodes increased.
	graph->at(ID)->push_back(edge);						// Insert the edge in the list of edges of ID
		numberOfEdges++;								// Increase the number of edges.
	updateReadLocations(edge);							// If the current edge contains some reads, then we need to update their location information.
	return true;
}


//=============================================================================
// Insert an edge in the graph.
//=============================================================================
bool OverlapGraph::insertEdge(Read *read1, Read *read2, UINT8 orient, UINT32 overlapOffset)
{
	Edge * edge1 = new Edge(read1,read2,orient,overlapOffset);			// Create a new edge in the graph to insert.
	// Set the overlap offset accordingly for the reverse edge. Note that read lengths are different.
	UINT32 overlapOffsetReverse = read2->getReadLength() + overlapOffset - read1->getReadLength();
	// If read lengths are the same. Then the reverse edge has the same overlap offset.
	Edge * edge2 = new Edge(read2,read1,twinEdgeOrientation(orient),overlapOffsetReverse);		// Create a new edge for the reverses string.

	edge1->setReverseEdge(edge2);		// Set the reverse edge pointer.
	edge2->setReverseEdge(edge1);		// Set the reverse edge pinter.
	insertEdge(edge1);					// Insert the edge in the overlap graph.
	insertEdge(edge2);					// Insert the edge in the overlap graph.
	return true;
}


//=============================================================================
// Orientation of a reverse edge;
// Twin edge of Orientation 0 = <-------< is Orientation 3 = >------->
// Twin edge of Orientation 1 = <-------> is Orientation 1 = <------->
// Twin edge of Orientation 2 = >-------< is Orientation 2 = >-------<
// Twin edge of Orientation 3 = >-------> is Orientation 0 = <-------<
//=============================================================================
UINT8 OverlapGraph::twinEdgeOrientation(UINT8 orientation)
{
	UINT8 returnValue;
	if(orientation == 0)
		returnValue = 3;
	else if(orientation == 1)
		returnValue = 1;
	else if(orientation == 2)
		returnValue = 2;
	else if(orientation == 3)
		returnValue = 0;
	else
		MYEXIT("Unsupported edge orientation.")
	return returnValue;
}


//=============================================================================
// Update location of reads for the new edge.
// The new edge many contain some reads. We need to update the location information of all such read.
//=============================================================================
bool OverlapGraph::updateReadLocations(Edge *edge)
{
	UINT64 distance = 0;
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++) 				// For each read in the edge
	{
		distance += edge->getListOfOverlapOffsets()->at(i);					// Distance of the read in the edge.
		Read *read = dataSet->getReadFromID(edge->getListOfReads()->at(i));	// Get the read object.
		if(edge->getListOfOrientations()->at(i) == 1)						// Orientation 1 means that the forward string of the read is contained in this edge.
		{
			read->getListOfEdgesForward()->push_back(edge);					// Insert the current edge in the list of forward edges.
			read->getLocationOnEdgeForward()->push_back(distance);			// Also insert the distance within the edge.
			read->getListOfEdgesForward()->resize(read->getListOfEdgesForward()->size());
			read->getLocationOnEdgeForward()->resize(read->getLocationOnEdgeForward()->size());
		}
		else																// Orientation 0 means that the reverser string of the read is contained in this edge.
		{
			read->getListOfEdgesReverse()->push_back(edge);					// Insert the current edge in the list of reverser edges.
			read->getLocationOnEdgeReverse()->push_back(distance);			// Also insert the distance withing the edge.
			read->getListOfEdgesReverse()->resize(read->getListOfEdgesReverse()->size());
			read->getLocationOnEdgeReverse()->resize(read->getLocationOnEdgeReverse()->size());
		}
	}
	return true;
}


//=============================================================================
// For each node in the graph, sort all its incident edges according to destination read ID.
//=============================================================================
void OverlapGraph::sortEdges()
{
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(!graph->at(i)->empty())
		{
			sort(graph->at(i)->begin(), graph->at(i)->end(), compareEdgeID);
		}
	}
}


//=============================================================================
// Contract composite paths in the overlap graph.
// u*-------->v>---------*w  => u*--------------------*w
// u*--------<v<---------*w  => u*--------------------*w
//=============================================================================
UINT64 OverlapGraph::contractCompositePaths(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	for(UINT64 index = 1 ; index < graph->size(); index++)
	{
		if(graph->at(index)->size() == 2) // Check if the node has only two edges.
		{
			Edge *edge1 = graph->at(index)->at(0);  // First edge.
			Edge *edge2 = graph->at(index)->at(1);  // Second Edge.
			// NEW: Will allow multiple edges between two nodes even before the flow.
			//if(flowComputed == true || !isEdgePresent(edge1->getDestinationRead()->getReadNumber(), edge2->getDestinationRead()->getReadNumber()))
				// Before flow is computed we do not insert multiple edges between the same endpoints.
				// CP: Before flow is computed (flowComputed == false), why checking if there is an edge between the two destination reads?
				// BH: Before the flow is computed we do not want to insert multiple edges between the same nodes. This condition is required for CS2 minimum Cost flow algorithm
				// CP: Why the destination reads, not the source reads??
				// BH: We do not want to remove loops (a,a).o you remove
				// CP: using the example above, don't you need to check if u and w are different reads or not?
				// BH: We do not need to check if u and w are different read or not.
			{
				if( matchEdgeType(edge1->getReverseEdge(), edge2) && edge1->getSourceRead() != edge1->getDestinationRead()) // One incoming edge and one outgoing edge.
				{
					mergeEdges(edge1->getReverseEdge(),edge2);	// Merge the edges.
					counter++;									// Counter how many edges merged.
				}
			}
		}

	}
	cout << setw(10) << counter << " composite Edges merged." << endl;
	CLOCKSTOP;
	return counter;
}


//=============================================================================
// Merge two edges in the overlap graph.
// CP: the flow of the new edge is the minimum flow of the two old edges and the flow of the new edge is deducted from those of the old edges
// CP: Remove the old edges that doesn't have any flow left, but keep the old edges that still have flow left.
//=============================================================================
bool OverlapGraph::mergeEdges(Edge *edge1, Edge *edge2)
{
	Edge *edgeForward = new Edge(); // New forward edge.
	Read *read1 = edge1->getSourceRead(), *read2 = edge2->getDestinationRead();
	Edge *edgeReverse = new Edge(); // New reverse edge.

	UINT8 orientationForward = mergedEdgeOrientation(edge1,edge2);			// Orientation of the forward edge.
	UINT8 orientationReverse = twinEdgeOrientation(orientationForward);		// Orientation of the reverse edge.

	vector<UINT64> * listReadsForward = new vector<UINT64>;					// List of reads in the forward edge.
	vector<UINT32> * listOverlapOffsetsForward= new vector<UINT32>;			// List of Overlaps in the forward edge.
	vector<UINT8> * listOrientationsForward = new vector<UINT8>;			// List of orientations in the forward edge.

	mergeList(edge1, edge2, listReadsForward, listOverlapOffsetsForward, listOrientationsForward); // Merge the lists from the two edges.
	// Make the forward edge
	edgeForward->makeEdge(read1,read2,orientationForward, edge1->getOverlapOffset() + edge2->getOverlapOffset(),
			listReadsForward, listOverlapOffsetsForward, listOrientationsForward);

	vector<UINT64> * listReadsReverse = new vector<UINT64>;					// List of reads in the reverse edge.
	vector<UINT32> * listOverlapOffsetsReverse= new vector<UINT32>;				// List of overlaps in the reverse edge.
	vector<UINT8> * listOrientationsReverse = new vector<UINT8>;			// List of orientations in the reverse edge.

	// Merge the lists from the two reverse edges.
	mergeList(edge2->getReverseEdge(),edge1->getReverseEdge(), listReadsReverse, listOverlapOffsetsReverse,listOrientationsReverse);
	// Make the reverse edge
	edgeReverse->makeEdge(read2, read1, orientationReverse, edge2->getReverseEdge()->getOverlapOffset() +
			edge1->getReverseEdge()->getOverlapOffset(), listReadsReverse, listOverlapOffsetsReverse, listOrientationsReverse);

	edgeForward->setReverseEdge(edgeReverse);		// Update the reverse edge pointer.
	edgeReverse->setReverseEdge(edgeForward);		// Update the reverse edge pointer.


	UINT32 flow = min(edge1->flow,edge2->flow);		// We take the minimum of the flow in the new edge.

	edgeForward->flow = flow;						// Modify the flow in the new forward edge.
	edgeReverse->flow = flow;						// Modify the flow in the new reverse edge.


	insertEdge(edgeForward);						// Insert the new forward edge in the graph.
	insertEdge(edgeReverse);						// Insert the new reverse edge in the graph.

	edge1->flow = edge1->flow - flow;				// Remove the used flow from edge1.
	edge1->getReverseEdge()->flow = edge1->flow;	// Remove the used flow from the reverse of edge1.

	edge2->flow = edge2->flow - flow;				// Remove the used flow from edge2.
	edge2->getReverseEdge()->flow = edge2->flow;	// Remove the used flow from teh reverse of edge2.

	if(edge1->flow == 0 || flow == 0)				// If no flow left in edge1
		removeEdge(edge1, false);					// edge1 is deleted from the graph.
	if(edge2->flow == 0 || flow == 0)				// If now flow left in edge2
		removeEdge(edge2, false);					// edge 2 is deleted from the graph.

	return true;

}


//=============================================================================
// Merge the list of reads, list of overlap offsets and list of orientations of two edges.
//=============================================================================
bool OverlapGraph::mergeList(const Edge *edge1, const Edge *edge2, vector<UINT64> *listReads, vector<UINT32> *listOverlapOffsets, vector<UINT8> *listOrientations)
{
	UINT64 sum = 0;
	for(UINT64 i = 0; i < edge1->getListOfOrientations()->size(); i++) 			// Take the list from edge1.
	{
		listReads->push_back(edge1->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge1->getListOfOverlapOffsets()->at(i));
		listOrientations->push_back(edge1->getListOfOrientations()->at(i));
		sum += edge1->getListOfOverlapOffsets()->at(i);
	}
	listReads->push_back(edge1->getDestinationRead()->getReadNumber()); 		// Insert the common node of the two edges

	listOverlapOffsets->push_back(edge1->getOverlapOffset() - sum);				// Get the overlap offset.

	if(edge1->getOrientation() == 1 || edge1->getOrientation() == 3)			// Orientation of the common node. Depends on the orientation of the edges.
		listOrientations->push_back(1);
	else
		listOrientations->push_back(0);
	for(UINT64 i = 0; i < edge2->getListOfOrientations()->size(); i++)			// take the list from edge2.
	{
		listReads->push_back(edge2->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge2->getListOfOverlapOffsets()->at(i));
		listOrientations->push_back(edge2->getListOfOrientations()->at(i));
	}
	return true;
}


//=============================================================================
// Orientation of the merged edge.
// e1(u,v) e2(v,w)
// Orientation e1 0 = u<-------<v		Orientation e2 0 = v<-------<w
// Orientation e1 1 = u<------->v		Orientation e2 1 = v<------->w
// Orientation e1 2 = u>-------<v		Orientation e2 2 = v>-------<w
// Orientation e1 3 = u>------->v		Orientation e2 3 = v>------->w
//
// 0 + 0 = u<-------<w
// 0 + 1 = u<------->w
//...................
//
//=============================================================================
UINT8 OverlapGraph::mergedEdgeOrientation(const Edge *edge1, const Edge *edge2)
{
	UINT8 or1 = edge1->getOrientation(), or2 = edge2->getOrientation(),returnValue;
	if(or1 == 0 && or2 == 0)
		returnValue = 0;
	else if(or1 == 0 && or2 == 1)
		returnValue = 1;
	else if(or1 == 1 && or2 == 2)
		returnValue = 0;
	else if(or1 == 1 && or2 == 3)
		returnValue = 1;
	else if(or1 == 2 && or2 == 0)
		returnValue = 2;
	else if(or1 == 2 && or2 == 1)
		returnValue = 3;
	else if(or1 == 3 && or2 == 2)
			returnValue = 2;
	else if(or1 == 3 && or2 == 3)
			returnValue = 3;
	else
	{
		cout<<(int)or1<<" "<<(int)or2<<endl;
		MYEXIT("Unable to merge.")
	}
	return returnValue;
}

//=============================================================================
// Remove an edge from the overlap graph.
//=============================================================================
bool OverlapGraph::removeEdge(Edge *edge, bool reportRemovedEdgeTag)
{
	// report edge to contig before removing it
	if (reportRemovedEdgeTag)
	{
		string s = getStringInEdge(edge); 	// get the string in the edge.
		(*contigSequences).push_back(s);
	}
	// remove edge
	removeReadLocations(edge);								// If the current edge contains some reads. We have to update their location formation.
	removeReadLocations(edge->getReverseEdge());			// Same for the reverse edge.
	Edge *twinEdge = edge->getReverseEdge();
	UINT64 ID1 = edge->getSourceRead()->getReadNumber(), ID2 = edge->getDestinationRead()->getReadNumber();  // Get the source and destation read IDs.
	for(UINT64 i = 0; i< graph->at(ID2)->size(); i++) // Delete the twin edge first.
	{
		if(graph->at(ID2)->at(i) == twinEdge)
		{
			delete graph->at(ID2)->at(i);
			graph->at(ID2)->at(i) = graph->at(ID2)->at(graph->at(ID2)->size()-1);
			graph->at(ID2)->pop_back();
			if(graph->at(ID2)->empty())
				numberOfNodes--;
			numberOfEdges--;
			break;
		}
	}
	for(UINT64 i = 0; i< graph->at(ID1)->size(); i++) // Delete the edge then.
	{
		if(graph->at(ID1)->at(i) == edge)
		{
			delete graph->at(ID1)->at(i);
			graph->at(ID1)->at(i) = graph->at(ID1)->at(graph->at(ID1)->size()-1);
			graph->at(ID1)->pop_back();
			if(graph->at(ID1)->empty())
				numberOfNodes--;
			numberOfEdges--;
			break;
		}
	}
	return true;
}


//=============================================================================
// Remove read mapping information of the current edge.
// This function is called when an edge is removed.
//=============================================================================
bool OverlapGraph::removeReadLocations(Edge *edge)
{
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++) 						// For each read in this edge.
	{
		Read *read = dataSet->getReadFromID(edge->getListOfReads()->at(i)); 		// Get the read object.
		for(UINT64 j = 0; j < read->getListOfEdgesForward()->size(); j++)  			// Check the list of edges that contain forward of this read.
		{
			if(read->getListOfEdgesForward()->at(j) == edge)						// Remove the current edge from the list.
			{
				read->getListOfEdgesForward()->at(j) = read->getListOfEdgesForward()->at(read->getListOfEdgesForward()->size()-1);
				read->getLocationOnEdgeForward()->at(j) = read->getLocationOnEdgeForward()->at(read->getLocationOnEdgeForward()->size()-1);

				read->getListOfEdgesForward()->pop_back();
				read->getLocationOnEdgeForward()->pop_back();

				read->getListOfEdgesForward()->resize(read->getListOfEdgesForward()->size());
				read->getLocationOnEdgeForward()->resize(read->getLocationOnEdgeForward()->size());
			}
		}

		for(UINT64 j = 0; j < read->getListOfEdgesReverse()->size(); j++) 			// Check the list of edges that contain reverse of this read.
		{
			if(read->getListOfEdgesReverse()->at(j) == edge) 						// Remove the current edge from the list.
			{
				read->getListOfEdgesReverse()->at(j) = read->getListOfEdgesReverse()->at(read->getListOfEdgesReverse()->size()-1);
				read->getLocationOnEdgeReverse()->at(j) = read->getLocationOnEdgeReverse()->at(read->getLocationOnEdgeReverse()->size()-1);

				read->getListOfEdgesReverse()->pop_back();
				read->getLocationOnEdgeReverse()->pop_back();

				read->getListOfEdgesReverse()->resize(read->getListOfEdgesReverse()->size());
				read->getLocationOnEdgeReverse()->resize(read->getLocationOnEdgeReverse()->size());
			}
		}
	}
	return true;
}



//=============================================================================
// Remove nodes with all simple edges and all same arrow type
//
// A node is all incoming edges or all outgoing edges is called a dead end node.
// While traversing the graph if we // enter such node, there is no way we can
// go out. So we remove such nodes from the graph. To remove the node all of
// its edges must be simple edges or very short edges (less than 10 read in it).
//=============================================================================
UINT64 OverlapGraph::removeDeadEndNodes(void)
{
	CLOCKSTART;
	vector <UINT64> listOfNodesToBeRemoved; // for saving nodes that should be deleted
	UINT64 edgesRemoved = 0;
	for(UINT64 i = 1; i < graph->size(); i++) // For each read.
	{
		if(!graph->at(i)->empty())	// If the read has some edges.
		{
			bool isDeadEnd = true;	// flag for dead end edge
			UINT64 inEdge = 0; 		// number of incoming edges to this node
			UINT64 outEdge = 0; 	// number of outgoing edges from this node

			for(UINT64 j=0; j < graph->at(i)->size(); j++) // For each edge
			{
				Edge * edge = graph->at(i)->at(j);
				// Break case:
				// 1. if composite edge is more than minReadsCountInEdgeToBeNotDeadEnd (deafult: 10)
				// 2. if composite edge is longer than minEdgeLengthToBeNotDeadEnd (default: 500)
				// 3. if the edge is loop for the current node
				// Then flag=1 and exit the loop

				// Case 1:
				if(edge->getListOfReads()->size() >= minReadsCountInEdgeToBeNotDeadEnd) {
					isDeadEnd = false;
					break;
				}
				// Case 2:
				if(edge->getEdgeLength() >= minEdgeLengthToBeNotDeadEnd) {
					isDeadEnd = false;
					break;
				}
				// Case 3:
				if(edge->getSourceRead()->getReadNumber() == edge->getDestinationRead()->getReadNumber())
				{
					isDeadEnd = false;
					break;
				}

				if(edge->getOrientation() == 0 || edge->getOrientation() == 1)
					inEdge++;			// the number of incoming edges to this node
				else
					outEdge++;			// the number of outgoing edges to this node
			}
			// if all edges are incoming or outgoing AND has less than 10 reads, this node is marked to be removed?
			// This node does not have any composite edge with more than 10 reads in it.
			if( isDeadEnd ) // If not break case
			{
				if( (inEdge > 0 && outEdge == 0) || (inEdge == 0 && outEdge > 0)) // only one type of simple edges
				{
					listOfNodesToBeRemoved.push_back(i);
				}
			}
		}
	}
	cout<< "Dead-end nodes detected: " << listOfNodesToBeRemoved.size() << endl;

    //  in below actually to remove the marked deadend nodes and all their edges.
    vector <Edge *> listOfEdges;
    for(UINT64 i = 0 ; i < listOfNodesToBeRemoved.size(); i++)
    {
        listOfEdges.clear();
        if(!graph->at(listOfNodesToBeRemoved.at(i))->empty())   // If the read has some edges.
        {
            edgesRemoved += graph->at(listOfNodesToBeRemoved.at(i))->size();
            for(UINT64 j=0; j < graph->at(listOfNodesToBeRemoved.at(i))->size(); j++) // For each edge
            {
                listOfEdges.push_back(graph->at(listOfNodesToBeRemoved.at(i))->at(j));
            }
            for(UINT64 j = 0; j< listOfEdges.size(); j++)
            {
                removeEdge(listOfEdges.at(j), reportRemovedEdgesToContigs);					// Remove all the edges of the current node.
            }
        }
    }
    cout<< "Dead-end nodes removed: " << listOfNodesToBeRemoved.size() << endl;
    cout<< "Dead-end edges removed: " << edgesRemoved << endl;

    CLOCKSTOP;
    return listOfNodesToBeRemoved.size();

}


//=============================================================================
// Simplify graph
//=============================================================================
bool OverlapGraph::simplifyGraph(void)
{
	UINT64 counter = 0;
	do
	{
		 counter  = removeDeadEndNodes();			// Remove dead-ends
		 counter += contractCompositePaths();		// Contract composite paths
		 counter += removeSimilarEdges();
		 // the two steps below requires flow to be computed
		 // if simplifyGraph is called in the unitig graph, these two functions will just return.
		 counter += reduceTrees();					// Remove trees.
		 counter += reduceLoops();					// Reduce loops

	} while (counter > 0);
	return true;
}


//=============================================================================
// Remove edges with similar endpoint in the overlap graph
// CP: Find pairs of edges with the same source and destination and representing similar sequences - 5% edit distance
// CP: Remove one of them and move its flow to the other edge
//=============================================================================
UINT64 OverlapGraph::removeSimilarEdges(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	// pairs of edges that are found, the two lists have the same length, a pair are in the same index
	vector <Edge *> listOfEdges1; // edges to be kept
	vector <Edge *> listOfEdges2; // edges to be removed

	vector <UINT64> listOfEditDistance;
	for(UINT64 i = 1; i < graph->size(); i++)	// For each node.
	{
		if(!graph->at(i)->empty())		// The node has some edges in the graph.
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)	// For all the edges.
			{
				Edge * e1 = graph->at(i)->at(j);
				UINT64 source1 = e1->getSourceRead()->getReadNumber(), destination1 = e1->getDestinationRead()->getReadNumber();
				if( source1 < destination1)
				{
					for(UINT64 k = j + 1; k < graph->at(i)->size(); k++)
					{
						Edge * e2 = graph->at(i)->at(k);
						UINT64 source2 = e2->getSourceRead()->getReadNumber(), destination2 = e2->getDestinationRead()->getReadNumber();
						if( source1 == source2 && destination1 == destination2)		// Check both edges starts and ends at same end points
						{
							if( abs((int)(e1->getOverlapOffset() - e2->getOverlapOffset())) <  (int)(e2->getOverlapOffset()/20) ) // The lengths are more than 95% similar
							{
								string s1 = getStringInEdge(e1);		// Get the string on the first edge.
								string s2 = getStringInEdge(e2);		// Get the string on the second edge.
								UINT64 editDistance = calculateEditDistance(s1,s2);		// Calculate the edit distance between the strings.
									// If the edit distance is less than 5% of the length of the shortest string.
								if( editDistance < min(e1->getOverlapOffset(), e2->getOverlapOffset())/20 )
								{
									UINT64 l;
									getBaseByBaseCoverage(e1);
									getBaseByBaseCoverage(e2);
									if(e1->coverageDepth < e2->coverageDepth ||			// first decide which edge to keep based on their coverage
																			// if the coverage are the same (most likely both coverage are zero), then decide based on the number of contained reads
									(e1->coverageDepth == e2->coverageDepth && e1->getListOfReads()->size() < e2->getListOfReads()->size()))
										// Will swap the edges based on coverage depth.
									{
										Edge *temp = e1;
										e1 = e2;
										e2 = temp;
									}
									for(l= 0; l <  listOfEdges1.size(); l++)	// Already in the list.
									{
										if(listOfEdges2.at(l) == e1 || listOfEdges2.at(l) ==  e2) // check if the edges are already used.
											break;
									}
									if(l ==  listOfEdges1.size())				// Not in the list. Add it in the list.
									{
										listOfEdges1.push_back(e1);					// We will keep this edge.
										listOfEdges2.push_back(e2);					// This edge will be deleted and the flow will be moved to the first edge.
										listOfEditDistance.push_back(editDistance);	// Also store the edit distance.
									}
								}
							}
						}
					}
				}
			}
		}
	}

	// Remove edges in listOfEdges2 and move their flow to their corresponding edge in listOfEdges1
	cout << listOfEdges1.size()<< " edges to remove" << endl;
	for(UINT64 i = 0; i < listOfEdges1.size(); i++)
	{
		cout << setw(10) << ++ counter << " removing edge ("<< setw(10) << listOfEdges1.at(i)->getSourceRead()->getReadNumber()<<","
				<< setw(10) << listOfEdges1.at(i)->getDestinationRead()->getReadNumber()<<") Lengths : " << setw(10)
				<< listOfEdges1.at(i)->getOverlapOffset() << " and " << setw(10) << listOfEdges2.at(i)->getOverlapOffset()
				<< " Flows: " << setw(3) << listOfEdges1.at(i)->flow << " and " << setw(3) << listOfEdges2.at(i)->flow
				<< " Edit Distance: " << setw(5) << listOfEditDistance.at(i) << " Reads: " << listOfEdges1.at(i)->getListOfReads()->size()
				<< " and " << listOfEdges2.at(i)->getListOfReads()->size() << endl;
		listOfEdges1.at(i)->flow += listOfEdges2.at(i)->flow;			// Move the flow of the delete edge to this edge.
		listOfEdges1.at(i)->getReverseEdge()->flow += listOfEdges2.at(i)->getReverseEdge()->flow;	// Same for the reverse edge.
		removeEdge(listOfEdges2.at(i), reportRemovedEdgesToContigs);
	}
	cout << counter << " edges removed." << endl;
	CLOCKSTOP;
	return counter;
}


//=============================================================================
// This function returns the string by overlapping the reads in an edge in the overlap graph
//=============================================================================
string OverlapGraph::getStringInEdge(const Edge *edge)
{
	// sequence of the source read
	string sourceRead = (edge->getOrientation() == 2 || edge->getOrientation() == 3)
			?  edge->getSourceRead()->getStringForward() : edge->getSourceRead()->getStringReverse();
	// sequence of the destination read
	string destinationRead = (edge->getOrientation() == 1 || edge->getOrientation() == 3)
			?  edge->getDestinationRead()->getStringForward() : edge->getDestinationRead()->getStringReverse();
	// sequence of the edge to be returned, starting with the source read
	string returnString = sourceRead;

	// sequence of the current edge in the middle of an edge
	string readTemp;
	// length of the previous read
	UINT64 previousLength = sourceRead.length();
	UINT64 substringLength = 0;

	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++)
	{
		readTemp = (edge->getListOfOrientations()->at(i) == 1)
				? dataSet->getReadFromID(edge->getListOfReads()->at(i))->getStringForward()
						: dataSet->getReadFromID(edge->getListOfReads()->at(i))->getStringReverse();
		// the length of the substring of the current read that is overhang or new from the previous read
		substringLength =  readTemp.length() + edge->getListOfOverlapOffsets()->at(i) - previousLength;
		// if the current read has no overlap with the previous read, they are connected by scaffolder and have a gap between them and insert an 'N' between the sequences
		if( edge->getListOfOverlapOffsets()->at(i) ==  previousLength)
			returnString = returnString + 'N' + readTemp.substr(readTemp.length() - substringLength, substringLength); // in this case, substringLength == readTemp.length()
		else
			returnString = returnString + readTemp.substr(readTemp.length() - substringLength, substringLength);
		previousLength = readTemp.length();


	}

	if(edge->getListOfReads()->empty()) // Simple edge
	{
		// if it's a simple edge, get the sequence from the source read and the destination read
		substringLength =  destinationRead.length() + edge->getOverlapOffset() - sourceRead.length();
		returnString = returnString + destinationRead.substr(destinationRead.length() - substringLength, substringLength);
	}
	else
	{
		// if it's a composite edge, add the string of the destination read
		substringLength = edge->getReverseEdge()->getListOfOverlapOffsets()->at(0);
		returnString = returnString + destinationRead.substr(destinationRead.length() - substringLength, substringLength);
	}

	return returnString;
}


//=============================================================================
// This function returns the edit distance between two strings.
// Code downloaded from http://rosettacode.org/wiki/Levenshtein_distance#C.2B.2B
//
// the Levenshtein distance is a metric for measuring the amount of difference between two sequences (i.e. an edit distance).
// The Levenshtein distance between two strings is defined as the minimum number of edits needed to transform one string into the other,
// with the allowable edit operations being insertion, deletion, or substitution of a single character.
// For example, the Levenshtein distance between "kitten" and "sitting" is 3, since the following three edits change one into the other,
// and there is no way to do it with fewer than three edits
//
//=============================================================================
UINT64 OverlapGraph::calculateEditDistance(const  string & s1, const string & s2)
{
	const UINT64 m(s1.size());
	const UINT64 n(s2.size());
	if( m==0 )
		return n;
	if( n==0 )
		return m;
	UINT64 *costs = new UINT64[n + 1];
	for( UINT64 k=0; k<=n; k++ )
		costs[k] = k;

	UINT64 i = 0;
	for ( std::string::const_iterator it1 = s1.begin(); it1 != s1.end(); ++it1, ++i )
	{
		costs[0] = i+1;
		UINT64 corner = i;
		UINT64 j = 0;
		for ( std::string::const_iterator it2 = s2.begin(); it2 != s2.end(); ++it2, ++j )
		{
			UINT64 upper = costs[j+1];
			if( *it1 == *it2 )
			{
				costs[j+1] = corner;
			}
			else
			{
				UINT64 t(upper<corner?upper:corner);
				costs[j+1] = (costs[j]<t?costs[j]:t)+1;
			}
			corner = upper;
		}
	}
	UINT64 result = costs[n];
	delete [] costs;
	//cout << s1 << endl << s2 << endl << result<< endl;
	return result;
}


//=============================================================================
// Calculate the coverage depth of an edge for every basepair and then update
// the Mean and SD of coverage depth in the edge. Only consider reads that are
// unique to the edge.
// CP: Calculate the variable, coverageDepth and SD, the Edge class
//=============================================================================
void OverlapGraph::getBaseByBaseCoverage(Edge *edge)
{
	vector<UINT64> * coverageBaseByBase = new vector<UINT64>;
	// Array lenght same as the string length in the edge.
	UINT64 length = edge->getOverlapOffset() + edge->getDestinationRead()->getReadLength();
	for(UINT64 i = 0; i <=length; i++)
	{
		coverageBaseByBase->push_back(0);				// At first all the bases are covered 0 times.
	}
	UINT64 overlapOffset = 0;
	// Increment the coverage of the section that each read covers,
	// NOT counting the source read and destination read because they are shared amony multiple edges
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++)	// For each read in the edge.
	{
		Read *read = dataSet->getReadFromID(edge->getListOfReads()->at(i));
		overlapOffset += edge->getListOfOverlapOffsets()->at(i);		// Where the current read starts in the string.
		UINT64 readLength = read->getReadLength();
		for(UINT64 j = overlapOffset; j < overlapOffset + readLength; j++)
		{
			coverageBaseByBase->at(j) += read->getFrequency();		// Increase the coverage of all bases by the frequency of the read.
		}
	}

	// clear the coverage of the section that is covered by non-unique reads that are contained by multiple edges
	overlapOffset = 0;
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++) // Scan the reads again.
	{
		Read *read = dataSet->getReadFromID(edge->getListOfReads()->at(i));
		overlapOffset += edge->getListOfOverlapOffsets()->at(i);
		UINT64 readLength = read->getReadLength();
		if(read->getListOfEdgesForward()->size() > 1)		// Clear the bases that are covered by reads apperaing in multiple places in the graph.
		{
			for(UINT64 j = overlapOffset; j < overlapOffset + readLength; j++)
			{
				coverageBaseByBase->at(j) = 0;
			}
		}
	}
	// clear the coverage of the section that is covered by the source read, because source read is present in multiple edges
	for(UINT64 i = 0; i < edge->getSourceRead()->getReadLength(); i++)
	{
		coverageBaseByBase->at(i) = 0;
	}

	// clear the coverage of the section that is covered by the destination read, because destination read is present in multiple edges
	for(UINT64 i = 0; i < edge->getDestinationRead()->getReadLength(); i++)	// Similarly clear the bases covered by the destination read.
	{
		coverageBaseByBase->at(coverageBaseByBase->size() - 1 - i) = 0;
	}

	UINT64 sum = 0, variance=0, count = 0, mean = 0, sd = 0;

	// calculate mean and standard deviation of coverage of all bases from non-zero counts
	for(UINT64 i = 0; i < coverageBaseByBase->size(); i++)
	{
		if(coverageBaseByBase->at(i) != 0)		// Count only the non-zero values.
		{
			sum += coverageBaseByBase->at(i);
			count++;
		}
	}
	if ( count != 0 )
	{
		mean = sum/count;		// Calculate the mean.

		for(UINT64 i = 0; i < coverageBaseByBase->size(); i++)
		{
			if(coverageBaseByBase->at(i) != 0)	// Calculate the variance.
			{
				variance += (mean - coverageBaseByBase->at(i)) * (mean - coverageBaseByBase->at(i));
			}
		}
		sd = sqrt(variance/count);	// Calculate the standard deviation.
	}
	// if no base of this edge is covered by unique reads, then the mean and SD are 0
	edge->coverageDepth = mean;		// Update the mean of the current edge.
	edge->SD = sd;					// Update the standard deviation of the current edge.
	delete coverageBaseByBase;
}


//=============================================================================
// This function removes in-trees and out-trees.
//=============================================================================
UINT64 OverlapGraph::reduceTrees(void)
{
	CLOCKSTART;
	if(this->flowComputed == false)
	{
		cout << "Flow not computed." << endl;
			return 0;
	}
	UINT64 NumOfInEdges, NumOfOutEdges, inFlow, outFlow, nodeMerged = 0;
	vector <Edge *> listOfInEdges, listOfOutEdges;
	for(UINT64 i = 0; i < graph->size(); i++)					// For each node in the graph
	{

		NumOfInEdges = 0; NumOfOutEdges = 0; inFlow = 0; outFlow = 0;
		listOfInEdges.clear(); listOfOutEdges.clear();
		for(UINT64 j = 0; j< graph->at(i)->size(); j++)
		{
				// Some conditions for not being considered as a tree
			if(graph->at(i)->at(j)->flow == 0 || graph->at(i)->at(j)->flow != graph->at(i)->at(j)->getReverseEdge()->flow
					|| graph->at(i)->at(j)->getSourceRead()->getReadNumber() == graph->at(i)->at(j)->getDestinationRead()->getReadNumber() )
					break;
			if(graph->at(i)->at(j)->getOrientation() == 0 || graph->at(i)->at(j)->getOrientation() == 1 )		// It is an in-edge
			{
				NumOfInEdges++;
				inFlow += graph->at(i)->at(j)->flow;															// Count the in flow
				listOfInEdges.push_back(graph->at(i)->at(j));
			}
			else																								// It is an out-edge
			{
				NumOfOutEdges++;
				outFlow += graph->at(i)->at(j)->flow;															// Count the out flow
				listOfOutEdges.push_back(graph->at(i)->at(j));
			}

			if(inFlow == outFlow && ( (NumOfInEdges == 1 && NumOfOutEdges > 1) || (NumOfInEdges > 1 && NumOfOutEdges == 1) ) )		// Either an in tree or an out tree
			{
				nodeMerged++;
				for(UINT64 k = 0; k < listOfInEdges.size(); k++)
				{
					for(UINT64 l = 0; l < listOfOutEdges.size(); l++)
					{
						mergeEdges(listOfInEdges.at(k)->getReverseEdge(), listOfOutEdges.at(l));
					}
				}
			}
		}
	}
	cout << setw(10) << nodeMerged << " trees removed." << endl;
	CLOCKSTOP;
	return nodeMerged;
}


//=============================================================================
// This function remove loops that can be traversed in only one way.
// a>--->b>--->b>--->c
//
// this function doesn't resolve loops that can be traversed in two ways:
// example: a<---<b<--->b>--->c
//=============================================================================
UINT64 OverlapGraph::reduceLoops(void)
{
	CLOCKSTART;
	if(this->flowComputed == false)
	{
		cout << "Flow not computed." << endl;
		return 0;
	}
	UINT64 counter = 0;
	Edge *ab,*bb,*bc;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(graph->at(i)->size() == 4) // only four edges. The loop is counted twice.
		{
			UINT64 loopCount = 0, incomingEdgeCount = 0, outgoingEdgeCount = 0;
			for(UINT64 j = 0; j< graph->at(i)->size(); j++)
			{
				if(graph->at(i)->at(j)->getDestinationRead()->getReadNumber() == i) // This is a loop
				{
					loopCount++;
					bb = graph->at(i)->at(j);
				}
				else if(graph->at(i)->at(j)->getOrientation() == 0 || graph->at(i)->at(j)->getOrientation() == 1) // incoming edge
				{
					incomingEdgeCount++;
					ab = graph->at(i)->at(j)->getReverseEdge();
				}
				else if(graph->at(i)->at(j)->getOrientation() == 2 || graph->at(i)->at(j)->getOrientation() == 3) // outgoing edge
				{
					outgoingEdgeCount++;
					bc = graph->at(i)->at(j);
				}
			}
			if(loopCount==2 && incomingEdgeCount == 1 && outgoingEdgeCount == 1)  // two in the loop and one incoming and one outgoing
			{
				cout<<"Loop found at node: " << i <<  " loop edge length: " << bb->getOverlapOffset()
						<< " flow: " << bb->flow << " Other edge lengths: " << ab->getOverlapOffset() << " and " << bc->getOverlapOffset() << endl;
				//cout << getStringInEdge(bb) << endl;
				if(bb->getOrientation() == 0)
				{
					counter++;
					mergeEdges(ab,bb->getReverseEdge());
				}
				else if(bb->getOrientation() == 3)
				{
					counter++;
					mergeEdges(ab,bb);
				}
				else			// Arrow in the edge is >---< or <---->. In this case it is not possible to decide which way to go.
				{
					cout << "Unable to reduce loop because of the edge type." << endl;
				}
			}
		}
	}
	cout <<" Loops removed: " << counter << endl; // Number of loop we were able to reduce
	CLOCKSTOP;
	return counter;
}


//=============================================================================
// Calculate min cost flow
//
// An sample input to the CS2 algorithm
//
// p min       3840      13449									// p min numberOfNodes numberOfEdges
// n          1         0										// initial flow of 0 in node 1, node 1 is the supersource
// n       3840         0										// initial flow of 0 in node 3840, which is the supersink (equal to the number of nodes)
// a       3840          1          1    1000000    1000000		// edge from supersink to supersource, LowBound(1), UpperBound(1000000), Cost per unit of flow(1000000)
// a          1          2          0    1000000          0		// edge from supersource to node 2, with the defined cost function
// this continues to
// connect each node to supersource and supersink
// connect every edge in the original graph
//=============================================================================
/*
bool OverlapGraph::calculateFlow(string inputFileName, string outputFileName)
{
	CLOCKSTART;
	for(UINT64 i = 1; i < graph->size(); i++) // clear all the flow
	{
		if(!graph->at(i)->empty())
		{
			for(UINT64 j = 0; j< graph->at(i)->size(); j++)
			{
				graph->at(i)->at(j)->flow =0;
			}

		}
	}
	// CP: the number of vertices for CS2 graph, two CS2 nodes for each read nodes plus supersource and supersink
	UINT64 V = numberOfNodes * 2 + 2;
	// CP: the number of edges for CS2 graph,
	// CP: 3 CS2 edges for each overlap graph edge, plus two edges for CS nodes to supersource and supersink, plus an edge between supersource and supersink
	UINT64 E = numberOfEdges * 3 + numberOfNodes * 4 + 1;
	UINT64 SUPERSOURCE = 1;
	UINT64 SUPERSINK = V;
	// For each edge in the overlap graph, we add 3 edges in the directed graph. Two nodes are created for each node in the original graph.
	// A super source and a super sink is added. Each node is connected to the super source and super sink.
	INT64 FLOWLB[3], FLOWUB[3], COST[3];			// Flow bounds and cost of the edges.
	// CP: comment on what's inputfile and output files.
	// CP: it's confusing to open outputFile with inputFileName
	// BH: This is an ouput file for this function, but an input file for CS2
	ofstream outputFile;
	outputFile.open(inputFileName.c_str());
	if(outputFile == NULL)
		MYEXIT("Unable to open file: "+inputFileName);

	stringstream ss;
	ss << "p min " << setw(10) << V << " " << setw(10) << E << endl;  	// Number of nodes and edges in the new graph.
	ss << "n " << setw(10) << SUPERSOURCE << setw(10) << " 0" << endl;	// Flow in the super source
	ss << "n " << setw(10) << SUPERSINK << setw(10) << " 0" << endl;	// Flow in the super sink.
	FLOWLB[0] = 1; FLOWUB[0] = 1000000; COST[0] = 1000000;
	// Add an edge from super sink to super source with very high cost.
	ss << "a " << setw(10) << SUPERSINK << " " << setw(10) << SUPERSOURCE << " "
			<< setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;

	// CP: this is a lookup table from the overlap graph node ID to the CS2 graph ID
	vector<UINT64> *listOfNodes = new vector<UINT64>;
	// CP: this is a lookup table from the CS2 graph ID to the overlap graph node ID
	vector<UINT64> *listOfNodesReverse = new vector<UINT64>;

	// For n nodes in the graph, CS2 requires that the nodes are numbered from 1 to n. In the overlap graph,
	// the nodes does not have sequencinal ID. We need to convert them to 1 - n

	// If the ID of a node in the original graph is 100 and directed graph is 5
	// Then listOfNodes->at(100) is equal to 5
	// and ListOfNodesReverse->at(5) is equal to 100.

	for(UINT64 i = 0; i <= graph->size(); i++)
	{
		listOfNodes->push_back(0);
		listOfNodesReverse->push_back(0);
	}

	// This loop set lower bound and upper bound of each node to super source and to super sink. All costs are 0.
	UINT64 currentIndex = 1;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(!graph->at(i)->empty()) // edges to and from the super source and super sink
		{
			listOfNodes->at(i) = currentIndex;					// Mapping between original node ID and cs2 node ID
			listOfNodesReverse->at(currentIndex) = i;			// Mapping between original node ID and cs2 node ID
			// CP: don't you create two CS2 node for each original node? where is the second CS2 node ID?
			// BH: Yes I created two nodes for CS2. For a node u in the listOfNodes. We created 2*u and 2*u+1 in for the cs2
			FLOWLB[0] = 0; FLOWUB[0] = 1000000; COST[0] = 0;
			ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << 2 * currentIndex << " " << setw(10) << FLOWLB[0]
			                  << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << 2 * currentIndex + 1 << " " << setw(10) << FLOWLB[0]
			                  << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			ss << "a " << setw(10) << 2 * currentIndex << " " << setw(10) << SUPERSINK << " " << setw(10) << FLOWLB[0]
			                  << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			ss << "a " << setw(10) << 2 * currentIndex + 1 << " " << setw(10) << SUPERSINK << " " << setw(10) << FLOWLB[0]
			                  << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			currentIndex++;
		}
	}

	// This loop converts the original bi-directed edges to directed edges (1 becomes 6).
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(!graph->at(i)->empty()) // edges to and from the super source and super sink
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)
			{

				// CP: u and v are the CS2 node IDs of the source node and destination node, respectively, of the edge
				Edge *edge = graph->at(i)->at(j);
				UINT64 u = listOfNodes->at(edge->getSourceRead()->getReadNumber());
				UINT64 v = listOfNodes->at(edge->getDestinationRead()->getReadNumber());

				// set the bound and cost here
				// if edge has more than 20 reads:
				//   FLOWLB[0] = 1; FLOWUB[0] = 1; COST[0] = 1;
				//   FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
				//   FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
				// else:
				//   FLOWLB[0] = 0; FLOWUB[0] = 1; COST[0] = 1;
				//   FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
				//   FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
				// cost function is set in such a way that for the composite edges with more that 20 reads in them we will push at least one units of flow.
				// Otherwise we set the lower bound of flow to 0, meaning that these edges might not have any flow at the endl.
				// The cost of pushing the first using of flow in very cheap and then we pay high price to push more than 1 units of flow.
				// This will ensure that we do not push flow were it is not necessary.
				// CP: if we need to change the cost function, we just need to change this function, right?
				// BH: Yes, we only need to change this function if we want to use different cost function.
				calculateBoundAndCost(edge, FLOWLB, FLOWUB, COST);


				if(u < v || (u == v && edge < edge->getReverseEdge()))
				{
					// Here for each edge we add three edges with different values of cost and bounds.
					// Total 6 edges considering the reverse edge too.
					// For details on how to convert the edges off different types please see my thesis.

					// BH: u1, u2, v1 and v2 are the actual CS2 node IDs
					UINT64 u1 = 2 * u, u2 = 2 * u + 1, v1 =  2 * v, v2 = 2 * v + 1;

					// Connect the edges of the original graph
					// for each orignal edge, we add six edges
					if(edge->getOrientation() == 0)
					{
						// first edge in the cost function, forward and reverse
						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[0] << " " << setw(10)
								<< FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[0] << " " << setw(10)
								<< FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						// second edge in the cost function, forward and reverse
						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[1] << " " << setw(10)
								<< FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[1] << " " << setw(10)
								<< FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						// third edge in the cost function, forward and reverse
						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[2] << " " << setw(10)
								<< FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[2] << " " << setw(10)
								<< FLOWUB[2] << " " << setw(10) << COST[2] << endl;
					}
					else if(edge->getOrientation() == 1)
					{
						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[0] << " " << setw(10)
								<< FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[0] << " " << setw(10)
								<< FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[1] << " " << setw(10)
								<< FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[1] << " " << setw(10)
								<< FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[2] << " " << setw(10)
								<< FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[2] << " " << setw(10)
								<< FLOWUB[2] << " " << setw(10) << COST[2] << endl;

					}
					else if(edge->getOrientation() == 2)
					{
						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[0] << " " << setw(10)
								<< FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[0] << " " << setw(10)
								<< FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[1] << " " << setw(10)
								<< FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[1] << " " << setw(10)
								<< FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[2] << " " << setw(10)
								<< FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[2] << " " << setw(10)
								<< FLOWUB[2] << " " << setw(10) << COST[2] << endl;

					}
					else if(edge->getOrientation() == 3)
					{
						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[0] << " " << setw(10)
								<< FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[0] << " " << setw(10)
								<< FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[1] << " " << setw(10)
								<< FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[1] << " " << setw(10)
								<< FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[2] << " " << setw(10)
								<< FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[2] << " " << setw(10)
								<< FLOWUB[2] << " " << setw(10) << COST[2] << endl;

					}
				}
			}
		}
	}
	outputFile << ss.str();		// Write the string in a file for CS2
	outputFile.close();

	// CP: clear ss, right?
	ss.str(std::string());

	// CP: what are  you doing here?
	// BH: CS2 requires the file name to be char * not string. We convert the string filename to char *
	char * inFile = new char[inputFileName.size() + 1];
	std::copy(inputFileName.begin(), inputFileName.end(), inFile);
	inFile[inputFileName.size()] = '\0';

	char * outFile = new char[outputFileName.size() + 1];
	std::copy(outputFileName.begin(), outputFileName.end(), outFile);
	outFile[outputFileName.size()] = '\0';


	cout << "Calling CS2" << endl;
	main_cs2(inFile,outFile);			// Call CS2
	cout << "CS2 finished" << endl;

	delete[] inFile;
	delete[] outFile;

	ifstream inputFile;
	inputFile.open(outputFileName.c_str());
	if(inputFile == NULL)
		MYEXIT("Unable to open file: "+outputFileName);

	UINT64 lineNum = 0;
	while(!inputFile.eof())
	{
		lineNum ++;
		UINT64 source, destination, flow;
		inputFile >> source >> destination >> flow;		// get the flow from CS2
		// CP: give an sample of the CS2 output
		//
		//From To Flow
		//1 2421 0 from node 1 to node 2421 with flow of 0
		//1 3 0	from node 1 to node 3 with flow of 0
		//
		if(source != SUPERSINK && source != SUPERSOURCE && destination != SUPERSOURCE && destination != SUPERSINK && flow!=0)
		{
			UINT64 mySource = listOfNodesReverse->at(source/2);				// Map the source to the original graph
			UINT64 myDestination = listOfNodesReverse->at(destination/2);	// Map the destination in the original graph
			Edge *edge = findEdge(mySource, myDestination);					// Find the edge in the original graph.
			edge->flow += flow;												// Add the flow in the original graph.
		}
	}
	inputFile.close();
	delete listOfNodes;
	delete listOfNodesReverse;
	this->flowComputed = true;
	CLOCKSTOP;
	return true;
}
*/


bool OverlapGraph::calculateFlowStream(void)
{
    CLOCKSTART;
    // Add super source and super sink nodes, add edge from super sink to super source with very big cost
    // Add edge from super source to every node in the graph, also from every node in the graph to the super sink
    // Every edge will be assigned a lower bound and an upper bound of the flow (capacity), and the cost associated with the edge
    // NEW CHANGES:
    // Now change the initial flow so that super source only connects to the source nodes (no in-edges only out-edges), and
    // only sink nodes connect to the super sink. This way the simple edges that are not really needed but that will make the paths longer will be kept.
    //
    UINT64 V = numberOfNodes + 2, E = numberOfEdges * 3 + numberOfNodes * 2 + 1 , SUPERSOURCE = 1, SUPERSINK = V;
    INT64 FLOWLB[3], FLOWUB[3], COST[3];            // Flow bounds and cost of the edges, cost function originally is a piecewise function with 3 segments
    stringstream ss;
    ss << "p min " << setw(10) << V << " " << setw(10) << E << endl;    // Problem description: Number of nodes and edges in the graph.
    ss << "n " << setw(10) << SUPERSOURCE << setw(10) << " 0" << endl;  // Flow in the super source
    ss << "n " << setw(10) << SUPERSINK << setw(10) << " 0" << endl;    // Flow in the super sink.

    // Flow lower bound and upper bound, and the cost for the first segment in the piecewise cost function
    FLOWLB[0] = 1; FLOWUB[0] = 1000000; COST[0] = 1000000;
    ss << "a " << setw(10) << SUPERSINK << " " << setw(10) << SUPERSOURCE << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl; // Add an edge from super sink to super source with very high cost (almost infinity), also at most can be used once


    // If the ID of a node in the original graph is 100 and directed graph is 5
    // Then listOfNodes->at(100) is equal to 5
    // and ListOfNodesReverse->at(5) is equal to 100.
    vector<UINT64> *listOfNodes = new vector<UINT64>;
    vector<UINT64> *listOfNodesReverse = new vector<UINT64>;

    for(UINT64 i = 0; i <= graph->size(); i++)      // For n nodes in the graph, CS2 requires that the nodes are numbered from 1 to n. In the overlap graph, the nodes does not have sequencinal ID. We need to convert them to 1 - n
    {
        listOfNodes->push_back(0);
        listOfNodesReverse->push_back(0);
    }

    // This loop set lower bound and upper bound from super source to each node, and from each node to super sink. All costs are 0.
    UINT64 currentIndex = 1;
    for(UINT64 i = 1; i < graph->size(); i++)
    {
    	// JJ added numInEdges at Read class.
        //if( (dataSet->getReadFromID(i)->numInEdges + dataSet->getReadFromID(i)->numOutEdges) != 0 ) // edges to and from the super source and super sink
    	// This is from the original omega
    	if(!graph->at(i)->empty()) // edges to and from the super source and super sink
    	{
            // FILE_LOG(logDEBUG4) << "Found node " << i << " corresponding to index " << currentIndex;
            listOfNodes->at(i) = currentIndex;                  // Mapping between original node ID and cs2 node ID
            listOfNodesReverse->at(currentIndex) = i;           // Mapping between original node ID and cs2 node ID
            FLOWLB[0] = 0; FLOWUB[0] = 1000000; COST[0] = 0;
            ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << currentIndex + 1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
            ss << "a " << setw(10) << currentIndex + 1 << " " << setw(10) << SUPERSINK << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
            currentIndex++;
        }
    }

    // This loop converts the original bi-directed edges to directed edges (1 becomes 6).
    // This loop set the lower and upper bounds of the flow in each edge, and the cost
    for(UINT64 i = 1; i < graph->size(); i++)
    {
        if(!graph->at(i)->empty()) // edges to and from the super source and super sink
        {
            for(UINT64 j = 0; j < graph->at(i)->size(); j++)
            {
                Edge *edge = graph->at(i)->at(j);
				UINT64 u = listOfNodes->at(edge->getSourceRead()->getReadNumber());
				UINT64 v = listOfNodes->at(edge->getDestinationRead()->getReadNumber());

                calculateBoundAndCost(edge, FLOWLB, FLOWUB, COST);  // Calculate bounds and cost depending on the edge property (number of reads in edge, or string length)

                // Here for each edge we add three edges with different values of cost and bounds, 3 pieces in the piecewise function
                ss << "a " << setw(10) << u + 1 << " " << setw(10) << v + 1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
                ss << "a " << setw(10) << u + 1 << " " << setw(10) << v + 1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
                ss << "a " << setw(10) << u + 1 << " " << setw(10) << v + 1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
            }
        }
    }

    stringstream oss;                       // stringstream to catch the output
    // FILE_LOG(logINFO) << "Calling CS2 for flow analysis";
    main_cs2(&ss, oss);         // Call CS2
    // FILE_LOG(logINFO) << "Flow analysis finished";


    string s, d, f;
    UINT64 lineNum = 0;
    while(!oss.eof())
    {
        lineNum ++;
        UINT64 source, destination, flow;
        oss >> source >> destination >> flow;       // get the flow from CS2

        if(source != SUPERSINK && source != SUPERSOURCE && destination != SUPERSOURCE && destination != SUPERSINK && flow!=0)
        {
            UINT64 mySource = listOfNodesReverse->at(source-1);             // Map the source to the original graph
            UINT64 myDestination = listOfNodesReverse->at(destination-1);   // Map the destination in the original graph
            Edge *edge = findEdge(mySource, myDestination); // Find the edge in the original graph.
            edge->flow += flow;                                             // Add the flow in the original graph.
            // FILE_LOG(logDEBUG2) << "Edge from " << mySource << " to " << myDestination << " has flow " << edge->flow;
        }
    }
    delete listOfNodes;
    delete listOfNodesReverse;
    this->flowComputed = true;
    CLOCKSTOP;
    return true;
}


//=============================================================================
// This function calculates the cost and bounds for an edge in the overlap graph.
// This function is very sensitive to the assembled contigs and to the cost function parameters
// BH: changing the bounds and threshold of number of nodes in the edge (here 20) many give use very wrong flow.
//
// CP: given an *edge, calculate and return FLOWLB, FLOWUB, and COST
// CP: FLOWLB, FLOWUB, and COST are all array of size 3 because each overlap graph edge is represented by 3 CS2 edges to define a cost function
//=============================================================================
bool OverlapGraph::calculateBoundAndCost(const Edge *edge, INT64* FLOWLB, INT64* FLOWUB, INT64* COST)
{
	for(UINT64 i = 0; i < 3; i++)		// For the simple edges we put very high cost
	{
		FLOWLB[i] = 0; FLOWUB[i] = 10; COST[i] = 500000;
	}

	if(!edge->getListOfReads()->empty()) // Composite Edge
	{
		// Case1: Composite Edge of at least minFlowReadCountThreshold (default: 20) reads. Must have at least one unit of flow.
		// Case2: Composite Edge length is longer than minFlowEdgeLengthThreshold (default: 1000)
		if(edge->getListOfReads()->size() >= minReadsCountInEdgeToBe1MinFlow || edge->getEdgeLength() > minEdgeLengthToBe1MinFlow)
		{
			// the first 1 flow must be used and has a cost of 1
			// each additional flow up to 8 flows has a cost of 100000

			// this edge carry the first 1 flow
			FLOWLB[0] = 1; FLOWUB[0] = 1; COST[0] = 1;
			// this edge carries the second 1 flow with high cost
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
			// this edge provides additional flow after the second flow
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
		}
		else // Short composite edge containing less than 20 reads. May have zero flow.
		{
			// the first 1 flow may not be required, but has a low cost of 1
			// each additional flow up to 8 flows has a cost of 100000

			// this edge carries the first 1 flow
			FLOWLB[0] = 0; FLOWUB[0] = 1; COST[0] = 1;
			// this edge carries the second unit of flow with high cost
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
			// this edge provides additional flow after the two units of flow.
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
		}
	}

	return true;
}


//=============================================================================
// Return edge between source and destination
//=============================================================================
Edge * OverlapGraph::findEdge(UINT64 source, UINT64 destination)
{
	for(UINT64 i = 0; i < graph->at(source)->size(); i++) // For the list of edges of the source node.
	{
		if(graph->at(source)->at(i)->getDestinationRead()->getReadNumber() == destination)	// check if there is an edge to destination
			return graph->at(source)->at(i);	// return the edge.
	}
	cout << "Check for error " << source << " to " << destination << endl;
	MYEXIT("Unable to find edge");
	return (graph->at(0)->at(0));
}


//=============================================================================
// Remove an all simple edge in the overlap graph that does not have any flow.
// CP: return the number of edges removed
//=============================================================================
UINT64 OverlapGraph::removeAllSimpleEdgesWithoutFlow()
{
	CLOCKSTART;
	vector <Edge *> listOfEdges;
	for(UINT64 i = 1; i < graph->size(); i++) // For each read.
	{
		if(!graph->at(i)->empty())	// If the read has some edges.
		{
			for(UINT64 j=0; j < graph->at(i)->size(); j++) // For each edge
			{
				Edge * edge = graph->at(i)->at(j);
			if(edge->getSourceRead()->getReadNumber() < edge->getDestinationRead()->getReadNumber()
					&& edge->getListOfReads()->empty()	// this condition check if this is a composite edge or a simply edge
					&& edge->flow == 0) // The edge is simple edge with no flow.
				{
					listOfEdges.push_back(edge); // Put in the list of edges to be removed.
				}
			}
		}
	}
	for(UINT64 i = 0 ; i < listOfEdges.size(); i++)
	{
		//removedEdges->push_back(listOfEdges.at(i));
		//(*removedEdges).push_back(*(listOfEdges.at(i)));
		//removeEdge(listOfEdges.at(i));		// remove the edges from the list.
		removeEdge(listOfEdges.at(i), reportRemovedEdgesToContigs);		// remove the edges from the list.
	}
	cout<<"Edges removed: " << listOfEdges.size() << endl;
	CLOCKSTOP;
	return listOfEdges.size();
}


/*=============================================================================
	This function prints the overlap graph in overlap_graph->gdl file. The graph can be viewed by
	aisee (free software available at http://www.aisee.com/)
	It also stores the contigs in a file.

graph: {
layoutalgorithm :forcedir
fdmax:704
tempmax:254
tempmin:0
temptreshold:3
tempscheme:3
tempfactor:1.08
randomfactor:100
gravity:0.0
repulsion:161
attraction:43
ignore_singles:yes
node.fontname:"helvB10"
edge.fontname:"helvB10"
node.shape:box
node.width:80
node.height:20
node.borderwidth:1
node.bordercolor:31
node: { title:"43" label: "43" }	// node, Title and label are both node ID 43 (and read number)
node: { title:"65" label: "65" }
............................................
............................................
// edges from source node 43 to destination node 32217, thickness of 3 means composite edge, thickness of 1 for simple edge
// edge type of backarrowstyle:solid arrowstyle:solid color: green is >----------------<
// edge type of arrowstyle:solid color: red is <----------------<
// edge type of arrowstyle: none color: blue  is <------------------->
// (1,0x,206,30) means (Flow, coverageDepth, OverlapOffset, numberOfReads)

edge: { source:"43" target:"32217" thickness: 3 backarrowstyle:solid arrowstyle:solid color: green label: "(1,0x,206,30)" }
edge: { source:"65" target:"38076" thickness: 3 arrowstyle:solid color: red label: "(0,0x,75,11)" }
edge: { source:"280" target:"47580" thickness: 3 arrowstyle: none color: blue label: "(0,0x,123,11)" }
}
=============================================================================*/


//=============================================================================
// Print graph
//=============================================================================
bool OverlapGraph::printGraph(string graphFileName, string contigFileName)
{
	CLOCKSTART;

	UINT64 thickness, highestDegree = 0, highestDegreeNode;
	ofstream graphFilePointer, contigFilePointer;
	graphFilePointer.open(graphFileName.c_str());
	contigFilePointer.open(contigFileName.c_str());

	// open files
	if(graphFilePointer == NULL)
		MYEXIT("Unable to open file: "+graphFileName);
	if(contigFilePointer == NULL)
		MYEXIT("Unable to open file: "+contigFileName);

	// graph file header
	graphFilePointer << "graph: {"
			<< endl <<  "layoutalgorithm :forcedir" << endl <<  "fdmax:704" << endl <<  "tempmax:254"
			<< endl <<  "tempmin:0" << endl <<  "temptreshold:3" << endl <<  "tempscheme:3" << endl <<  "tempfactor:1.08"
			<< endl <<  "randomfactor:100" << endl <<  "gravity:0.0" << endl <<  "repulsion:161" << endl <<  "attraction:43"
			<< endl <<  "ignore_singles:yes" << endl <<  "node.fontname:\"helvB10\"" << endl << "edge.fontname:\"helvB10\""
			<< endl <<  "node.shape:box" << endl <<  "node.width:80" << endl <<  "node.height:20" << endl <<  "node.borderwidth:1"
			<< endl <<  "node.bordercolor:31" << endl;
	for(UINT64 i = 1; i<= dataSet->getNumberOfUniqueReads(); i++)
		if(!graph->at(i)->empty())
			graphFilePointer << "node: { title:\""<<i<<"\" label: \"" << i << "\" }" << endl;

	// graph file
	for(UINT64 i = 1; i<= dataSet->getNumberOfUniqueReads(); i++)
	{
		if(!graph->at(i)->empty())
		{
			if(graph->at(i)->size() > highestDegree)
			{
				highestDegree = graph->at(i)->size();
				highestDegreeNode = i;
			}

			for(UINT64 j=0; j < graph->at(i)->size(); j++)
			{
				Edge * e = graph->at(i)->at(j);
				UINT64 source = e->getSourceRead()->getReadNumber(), destination = e->getDestinationRead()->getReadNumber();
				if(source < destination || (source == destination && e < e->getReverseEdge()) )
				{
					string s = getStringInEdge(e); 	// get the string in the edge.
					(*contigSequences).push_back(s);

					thickness = e->getListOfReads()->empty() ? 1: 3;
					if(e->getOrientation() == 0)
						graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination
						<< "\" thickness: " << thickness << " arrowstyle: none backarrowstyle: solid color: red label: \"("
						<< e->flow << "," << e->coverageDepth << "x," << e->getOverlapOffset() << "," << e->getListOfReads()->size()
						<< ")\" }" << endl;
					else if(e->getOrientation() == 1)
						graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: "
						<< thickness << " backarrowstyle:solid arrowstyle:solid color: green label: \"(" << e->flow << ","
						<< e->coverageDepth <<  "x," << e->getOverlapOffset() << "," << e->getListOfReads()->size()
						<< ")\" }" << endl;
					else if(e->getOrientation() == 2)
						graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: "
						<< thickness << " arrowstyle: none color: blue label: \"(" << e->flow << "," << e->coverageDepth << "x,"
						<< e->getOverlapOffset() << "," << e->getListOfReads()->size()
						<< ")\" }" << endl;
					else if(e->getOrientation() == 3)
						graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: "
						<< thickness << " arrowstyle:solid color: red label: \"(" << e->flow << "," << e->coverageDepth <<  "x,"
						<< e->getOverlapOffset() << "," << e->getListOfReads()->size()
						<< ")\" }" << endl;
				}
			}
		}
		else if (reportUnoverlappedReadsToContigs)
		{
			Read *read;
			read = dataSet->getReadFromID(i);
			string s = read->getStringForward();
			(*contigSequences).push_back(s);
		}
	}
	graphFilePointer << "}";
	cout << "Aisee graph written." << endl;

	sort((*contigSequences).begin(),(*contigSequences).end(),byStringLength); // Sort the contigs by their length.
	reverse((*contigSequences).begin(), (*contigSequences).end());

	// report contig
	UINT64 sum = 0;
	UINT64 contigId = 1;
	for(UINT64 i = 0; i < (*contigSequences).size(); i++) 		// Store the contigs in a file.
	{
		string s = (*contigSequences)[i];

		if (s.size() >= minReadsLengthTobeReported)			// only report contigs length >= minReadsLengthTobeReported
		{
			contigFilePointer << ">"<< contigId << endl;
			sum += s.length();
			UINT64 start=0;
			do
			{
				contigFilePointer << s.substr(start,100) << endl;  // save 100 BP in each line.
				start+=100;
			} while (start < s.length());
			contigId ++;
		}
	}
	cout << "Print contig 1 is done!" << endl;

	// Print some statistics about the graph.
	cout<< "Total contig length: " << sum <<" BP" << endl;
	cout<< "Number of Nodes in the graph: " << getNumberOfNodes() << endl;
	cout<< "Number of Edges in the graph: " << getNumberOfEdges()/2 << endl;

	UINT64 simEdges = 0, comEdges = 0, inEdges = 0, outEdges = 0;
	for(UINT64 i=0; i < graph->at(highestDegreeNode)->size(); i++)
	{
		if(graph->at(highestDegreeNode)->at(i)->getListOfReads()->empty())
			simEdges++;
		else
			comEdges++;
		if(graph->at(highestDegreeNode)->at(i)->getOrientation() == 0 || graph->at(highestDegreeNode)->at(i)->getOrientation() == 1)
			inEdges++;
		else
			outEdges++;
	}
	// Print some more statistics on the node with highest degree.
	cout<< "Highest Degree Read " << highestDegreeNode << " has " << highestDegree << " neighbors." << endl;
	cout << "In Edges: " << inEdges << " Out Edges: " << outEdges << " Simple Edges: " << simEdges << " Composite Edges: " << comEdges << endl;
	cout<< "String: " << dataSet->getReadFromID(highestDegreeNode)->getStringForward() << endl;

	graphFilePointer.close();
	contigFilePointer.close();

	CLOCKSTOP;
	return true;
}
