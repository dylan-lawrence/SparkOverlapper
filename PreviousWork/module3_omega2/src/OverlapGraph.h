#ifndef OVERLAPGRAPH_H
#define OVERLAPGRAPH_H

//============================================================================
// Name        : OverlapGraph.h
// Author      : Tae-Hyuk (Ted) Ahn
// Version     : V0.0.1
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : OverlapGraph header file
//============================================================================

#include "Config.h"
#include "DataSet.h"
#include "Edge.h"

class OverlapGraph
{
	private:
		DataSet *dataSet; 							// Pointer to the dataset containing all the reads. this is NOT modified here

        /*
         *  the overlap graph is a vector of nodes (i.e. reads) with their vectors of incident edges
         *  the readNumber of a read is equal to the index of the node in this vector. Because readNumber starts from 1, the
         *  to iterate through all edges
         *      for(UINT64 i = 1; i < graph->size(); i++) // For each read.
         *          for(UINT64 j=0; j < graph->at(i)->size(); j++) // For each edge
         *              Edge * edge = graph->at(i)->at(j);
         */
		vector< vector<Edge *> * > *graph;			// Adjacency list of the graph.
		vector <string> *contigSequences;			// List of removed edges.

		UINT64 numberOfNodes;						// Number of nodes in the overlap graph.
		UINT64 numberOfEdges;						// Number of edges in the overlap graph.
		bool flowComputed;							// Flag to check wheather the flow is computed or not.

		// Insert an edge in the overlap graph
		bool insertEdge(Edge * edge);

		// Insert an edge in the overlap graph
		bool insertEdge(Read *read1, Read *read2, UINT8 orient, UINT32 overlapOffset);

		// Orientation of the reverse edge.
		UINT8 twinEdgeOrientation(UINT8 orientation);

		// Merge two edges in the  overlap graph.
		bool mergeEdges(Edge *edge1, Edge *edge2);

		// When we merge two edges, we need to merge the list of ordered reads, their overlap offsets and orientations.
		bool mergeList(const Edge *edge1, const Edge *edge2, vector<UINT64> *listReads, vector<UINT32> *listOverlapOffsets, vector<UINT8> * ListOrientations);

		// Orientation of the edge when two edges are merged.
		UINT8 mergedEdgeOrientation(const Edge *edge1, const Edge *edge2);

		// Remove an edge from the overlap graph.
		bool removeEdge(Edge *edge, bool reportRemovedEdgeTag);

		// Remove the location of all the reads from the current edge. This function is called when an edge is removed.
		bool removeReadLocations(Edge *edge);

		// Update the location of all the reads in the current edge. This function is called when a new edge is inserted.
		bool updateReadLocations(Edge *edge);

		// Contract composite paths in the overlap graph.
		UINT64 contractCompositePaths(void);

		// Remove dead-ends from the overlap graph.
		UINT64 removeDeadEndNodes(void);

		// remove multi-edges with similar strings
		UINT64 removeSimilarEdges(void);

		// Get the string in an edge by overlapping the ordered reads in the edge.
		string getStringInEdge(const Edge *edge);

		// Find the edit distance between two strings. Used by removeSimilarEdges()
		UINT64 calculateEditDistance(const string & s1, const string & s2);

		// Get the coverage Mean and SD of an edge. Only considering the unique reads.
		void getBaseByBaseCoverage(Edge *edge);

		// Remove trees in the overlap graph.
		UINT64 reduceTrees(void);

		// Loops that can be traversed only one way
		UINT64 reduceLoops(void);

		// Calculate bounds and costs of flow for minimum cost flow in the overlap graph.
		bool calculateBoundAndCost(const Edge *edge, INT64* FLOWLB, INT64* FLOWUB, INT64* COST);

		// Find an edge from source to destination in the overlap graph.
		Edge *findEdge(UINT64 source, UINT64 destination);

	public:

		// Default constructor
		OverlapGraph(void);

		// Destructor
		~OverlapGraph();

		UINT64 getNumberOfEdges(void){return numberOfEdges;}		// Get the number of edges in the overlap graph.
		UINT64 getNumberOfNodes(void){return numberOfNodes;}		// Get the number of nodes in the overlap graph.

		// Set the dataset pointer.
		bool setDataSet(DataSet *dataset){dataSet=dataset;  return true;}

		// Build the overlap graph from read/edge files.
		bool buildOverlapGraphFromFiles(const vector<string> &readFilenameList, const vector<string> &edgeFilenameList);

		// Load read file
		//UINT64 loadReadFile(const string &readFilename, vector<UINT64> &readIDListInFile, UINT64 &readID);	// return count (the number of reads in the file)
		UINT64 loadReadFile(const string &readFilename, unordered_map<UINT64, UINT64> &readIDMap, UINT64 &readID);	// return count (the number of reads in the file)

		// Load edge file
		//void loadEdgeFile(const string &readFilename, vector<UINT64> &readIDListInFile);
		void loadEdgeFile(const string &readFilename, unordered_map<UINT64, UINT64> &readIDMap);

		// Sort edges of each read based on ID of the destination read. this is only for ordering edges for convenience in the output file
		void sortEdges();

		// Some simple simplification.
		bool simplifyGraph(void);

		// Calculate the minimum cost flow of the overlap graph using file
		// bool calculateFlow(string inputFileName, string outputFileName);
		bool calculateFlowStream(void);

		// Remove all simple edges without flow
		UINT64 removeAllSimpleEdgesWithoutFlow();

		// Store the overlap graph for visual display and also store the contigs/scaffods in a file.
		bool printGraph(string graphFileName, string contigFileName);

};



#endif /* OVERLAPGRAPH_H */
