//============================================================================
// Name        : main.cpp
// Author      : Tae-Hyuk (Ted) Ahn
// Version     : V0.0.1
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Main code
//============================================================================

#include "Config.h"
#include "DataSet.h"
#include "OverlapGraph.h"


int main(int argc, char **argv) {
	CLOCKSTART;

    // Parse command line options:
    if(!Config::setConfig(argc, argv)){
        cout << "Error: wrong configurations" << endl;
        return false;
    }

    // Get read/edge file(s) list
    vector<string> readFilenameList = Config::getReadFilenames();
    vector<string> edgeFilenameList = Config::getEdgeFilenames();

    // Get output filename prefix
    string outputFilenamePrefix = Config::getOutputFilenamePrefix();

    // Initialize
    DataSet *dataSet = new DataSet();
    OverlapGraph *overlapGraph = new OverlapGraph();

    // Set the dataset pointer.
    overlapGraph->setDataSet(dataSet);

    // Build overlap graph from read/edge files(s)
    overlapGraph->buildOverlapGraphFromFiles(readFilenameList, edgeFilenameList);
    overlapGraph->sortEdges();

    // Simplify the graph before the flow.
    overlapGraph->simplifyGraph();

	// Flow analysis
    // overlapGraph->calculateFlow(outputFilenamePrefix+"_flow.input", outputFilenamePrefix+"_flow.output");
    overlapGraph->calculateFlowStream();
	overlapGraph->removeAllSimpleEdgesWithoutFlow();

    // Simplify the graph before the flow.
    overlapGraph->simplifyGraph();

	// Print graph format file and contigs
	overlapGraph->printGraph(outputFilenamePrefix+"_graph1.gdl", outputFilenamePrefix+"_contigs1.fasta");

    CLOCKSTOP;
	return 0;
}
