//============================================================================
// Name        : Config.cpp
// Author      : Tae-Hyuk (Ted) Ahn
// Version     : V0.0.1
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Config cpp file
//============================================================================

#include "Config.h"


//=============================================================================
// Default constructor
//=============================================================================
Config::Config() {
}


//=============================================================================
// Default destructor
//=============================================================================
Config::~Config() {
}

//=============================================================================
//Here is to create and initialize the static members
//=============================================================================
vector<string> Config::readFilenamesList;
vector<string> Config::edgeFilenamesList;
string Config::outFilenamePrefix = "new_omega_out";

// global variables with default value
unsigned int minReadsCountInEdgeToBeNotDeadEnd = 10;
unsigned int minEdgeLengthToBeNotDeadEnd = 500;
unsigned int minReadsCountInEdgeToBe1MinFlow = 20;
unsigned int minEdgeLengthToBe1MinFlow = 1000;
unsigned int minReadsLengthTobeReported = 200;
bool reportRemovedEdgesToContigs = false;
bool reportUnoverlappedReadsToContigs = false;


//=============================================================================
// print help usage
//=============================================================================
void Config::printHelp()
{
    cout << endl
         << "  Usage:" << endl
         << "    new_omega [OPTION]...<PARAM>..." << endl
         << endl
         << "  <PARAM>" << std::endl
         << "    -f\t contained read reduction read filename(s) (comma separated fasta/fastq)" << endl
         << "    -e\t ovelapped edge property graph filename(s) (comma separated edge list)" << endl
         << endl
         << "  [OPTION]" << std::endl
         << "    -h/--help\t only print out the help contents" << endl
         << "    -o\t all output filename prefix (default: new_omega_out)" << endl
         << "    -mcd\t Minimum reads Count in an edge regarding to be not a Dead-end (default: 10)" << endl
         << "    -mld\t Minimum edge Length regarding to be not a Dead-end (default: 500)" << endl
         << "    -mcf\t Minimum reads Count in an edge regarding to be 1 minimum Flow (default: 20)" << endl
         << "    -mlf\t Minimum edge Length regarding to be 1 minimum Flow (default: 1000)" << endl
         << "    -mlr\t Minimum reads Length to be Reported (default: 200)" << endl
         << "    -rre\t Report Removed Edges to contigs (default: false)" << endl
         << "    -rur\t Report Unoverlapped Reads to contigs (default: false)" << endl
         << endl;
}

bool Config::setConfig(int argc, char **argv)
{
	vector<string> argumentsList;
	cout << "PRINTING ARGUMENTS" << endl;
	for(int i = 0; i < argc; i++)
	{
		cout << argv[i] << ' ';
	}
	cout << endl;
	while(argc--)
			argumentsList.push_back(*argv++);

	if(argumentsList.size() == 1)
	{
        Config::printHelp();
		return false;
	}

	for(UINT64 i = 1; i <= argumentsList.size()-1; i++)
	{
		// -h/--help: only print out the help contents
		if (argumentsList[i] == "-h" || argumentsList[i] == "--help")
	    {
	        Config::printHelp();
			exit(0);
	    }
		// -f: contained read reduction read filename(s) (comma separated fasta/fastq)
		else if (argumentsList[i] == "-f")
	    {
            string inputFilenames=argumentsList[++i];
            stringstream ss(inputFilenames);
            string item;
    
            while (getline(ss, item, ','))
            {
                Config::readFilenamesList.push_back(item);
            }
	    }
		// -e: ovelapped edge property graph filename(s) (comma separated edge list)
		else if (argumentsList[i] == "-e")
	    {
            string inputFilenames=argumentsList[++i];
            stringstream ss(inputFilenames);
            string item;
    
            while (getline(ss, item, ','))
            {
                Config::edgeFilenamesList.push_back(item);
            }
	    }
		// -o: all output filename prefix (default: new_omega_out)
		else if (argumentsList[i] == "-o")
        {
			Config::outFilenamePrefix = argumentsList[++i];
        }
        // -mcd: Minimum reads Count in an edge regarding to be not a Dead-end (default: 10)
		else if (argumentsList[i] == "-mcd")
        {
			istringstream iss_a(argumentsList[++i]);
			iss_a >> minReadsCountInEdgeToBeNotDeadEnd;
        }
        // -mld: Minimum edge Length regarding to be not a Dead-end (default: 500)
		else if (argumentsList[i] == "-mld")
        {
			istringstream iss_b(argumentsList[++i]);
			iss_b >> minEdgeLengthToBeNotDeadEnd;
        }
        // -mcf: Minimum reads Count in an edge regarding to be 1 minimum Flow (default: 20)
		else if (argumentsList[i] == "-mcf")
        {
			istringstream iss_c(argumentsList[++i]);
			iss_c >> minReadsCountInEdgeToBe1MinFlow;
        }
        // -mlf: Minimum edge Length regarding to be 1 minimum Flow (default: 1000)
		else if (argumentsList[i] == "-mlf")
        {
			istringstream iss_d(argumentsList[++i]);
			iss_d >> minEdgeLengthToBe1MinFlow;
        }
        // -mlr: Minimum reads Length to be Reported (default: 200)
		else if (argumentsList[i] == "-mlr")
        {
			istringstream iss_d(argumentsList[++i]);
			iss_d >> minReadsLengthTobeReported;
        }
        // -rre\t Report Removed Edges to contigs (default: false)
		else if (argumentsList[i] == "-rre")
        {
			reportRemovedEdgesToContigs = true;
        }
        // -rur\t Report Unoverlapped Reads to contigs (default: false)
		else if (argumentsList[i] == "-rur")
        {
			reportUnoverlappedReadsToContigs = true;
        }
		else
        {
            Config::printHelp();
            return false;
        }
    }
	return true;
}


//=============================================================================
// Get read filenames
//=============================================================================
vector<string> Config::getReadFilenames()
{
    return Config::readFilenamesList;
}


//=============================================================================
// Get edge filenames
//=============================================================================
vector<string> Config::getEdgeFilenames()
{
    return Config::edgeFilenamesList;
}


//=============================================================================
// Get output prefix name
//=============================================================================
string Config::getOutputFilenamePrefix()
{
   return Config::outFilenamePrefix;
}

