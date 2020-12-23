// MultipleMapping.cpp : Defines the entry point for the application.
//

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <iterator>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <zlib.h>

#include "BITOps.h"
#include "LengthParse.h"
#include "MMCommon.h"
#include "ReadParse.h"
#include "BedgraphWrite.h"
#include "ReweightIterator.h"

using namespace std;
using namespace std::chrono;

//Define global variables

outStream outlog = { NULL };
treesVec tree1, tree1neg;
treesDoub tree2, tree2neg;
vector< vector<readMap> > reads_vector;
map<string, int> chrom_to_counter;
vector<string> counter_to_chrom, reads_str_vector;
vector<int> counter_to_length;
bool stranded, readoutput;
float fitrate;
ogzstream read_out_unit;

//Off-the-shelf solution to check whether the filename ends with the extension

bool hasEnding(const string& filename, const string& extension)
{
	if (filename.length() >= extension.length())
	{
		return (0 == filename.compare(filename.length() - extension.length(), extension.length(), extension));
	}
	else
	{
		return false;
	}
}

//Manages the procedure from reading the read files, running the iterator program, and bedgraph output, adjusting if this is part of a cross-validation.

void readThroughBedgraph(const vector<string> rfils, int iterations, float fixation, string output_prefix, float score_min, int crossval, int cval, bool contout, bool cross, bool onsa, int maxaligns)
{
	string read_ids = "";

	//id_counter serves to keep track of the available index of the reads vector so we don't have to keep calling size; initially set to 0 in the parser.
	int id_counter = -1;
	//raw_counter keeps track of the number of reads that we have parsed through, mostly for the purposes of cross-validation apportionment.
	int raw_counter = 0;

	for (unsigned int fils = 0; fils < rfils.size(); fils++)
	{
		outlog << "Now parsing read file " << rfils[fils] << "\n";
		auto start = high_resolution_clock::now();

		//Check if the file ends with ".gz" i.e. whether it is a compressed file or not.
		bool checkgz = hasEnding(rfils[fils], ".gz");

		if (checkgz)
		{
			//Using igzstream and appropriate function if compressed file.
			igzstream reads_file(rfils[fils].c_str());
			if (!reads_file) { cerr << "File " << rfils[fils] << " could not be opened for reading."; return; }
			if (stranded) { parseReadsFileStranded(reads_file, read_ids, id_counter, raw_counter, crossval, cval, maxaligns, score_min, onsa); }
			else { parseReadsFile(reads_file, read_ids, id_counter, raw_counter, crossval, cval, maxaligns, score_min, onsa); }
			reads_file.close();
		}
		else
		{
			//Using ifstream and appropriate function if not compressed file.
			ifstream reads_file(rfils[fils]);
			if (!reads_file) { cerr << "File " << rfils[fils] << " could not be opened for reading."; return; }
			if (stranded) { parseReadsFileStranded(reads_file, read_ids, id_counter, raw_counter, crossval, cval, maxaligns, score_min, onsa); }
			else { parseReadsFile(reads_file, read_ids, id_counter, raw_counter, crossval, cval, maxaligns, score_min, onsa); }
			reads_file.close();
		}

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<seconds>(stop - start);
		outlog << "Finished parsing read file " << rfils[fils] << " in " << duration.count() << " seconds\n\n";
	}
	read_ids.clear();

	//reweight_prefix is going to be used for bedgraph files, including continuous output and/or final output.
	string reweight_prefix = output_prefix;
	if (cross)
	{
		if (onsa) { reweight_prefix = reweight_prefix + "_" + to_string(crossval) + "-fold-only-" + to_string(cval+1); }
		else { reweight_prefix = reweight_prefix + "_" + to_string(crossval) + "-fold-sans-" + to_string(cval+1); }
	}

	if (contout)
	{
		outlog << "Writing bedgraph for iteration 0\n\n";
		string filename = reweight_prefix + "_iteration-0";
		writeBedgraphOutput(filename);
	}

	//Runs the reweight iterator function
	reweightIterator(iterations, fixation, reweight_prefix, contout);

	//Assesses whether a continous output iteration would have already output the final bedgraph, in which case final writing can be skipped
	bool finalwrite = iterations < 10 || iterations % 10 == 0;

	if (!contout || !finalwrite)
	{
		writeBedgraphOutput(reweight_prefix);
	}
}

int main(int argc, char* argv[])
{
	int opt;
	int iterations = 1, maxaligns = 50, crossval = 1;
	float fixation = 0, score_min = 0;
	string length_name, output_prefix;
	bool contout = false;
	stranded = false;
	readoutput = false;
	fitrate = 1;

	string helpmessage = "MultiMap for analysis of ambiguously mapping reads in ChIP-seq. Version 0.2.0.\n\n"
			"Usage: MultiMap [options] [bed or bed.gz file input(s)]\n\n"
			"Required options:\n"
			"-g : Genome length file listing all chromosomes and lengths. Chromosomes will appear in this order in output file.\n"
			"-o : Output prefix used for output bedgraph and log files.\n\n"
			"Non-required options:\n"
			"-i : Number of iterations. Default 1.\n"
			"-m : Maximum number of alignments for a read to be processed. Default 50.\n"
			"-s : Minimum score for Bowtie2 display. Default 0 (unscored).\n"
			"-v : N for N-fold cross-validation. Default 1 (no cross-validation).\n"
			"-c : Flag for continuous output bedgraphs. Default off.\n"
			"-S : Flag for strand-specific mode. Default off.\n"
			"-r : Flag for read output mode with weights. Default off.\n"
			"-l : Rate of fitting in reweighting. Default 1.\n"
			"-h : Display help message.\n\n"
			"Developed by Rohan Shah (rohanshah@uchicago.edu).\n";

	//Parse command options

	while ((opt = getopt(argc, argv, "g:o:i:x:m:s:v:l:cSrh")) != -1)
	{
		switch (opt)
		{
		case 'g':
			length_name = optarg;
			break;
		case 'o':
			output_prefix = optarg;
			break;
		case 'i':
			iterations = stoi(optarg);
			break;
		case 'x':
			fixation = stof(optarg);
			break;
		case 'm':
			maxaligns = stoi(optarg);
			break;
		case 's':
			score_min = stof(optarg);
			break;
		case 'v':
			crossval = stoi(optarg);
			break;
		case 'l':
			fitrate = stof(optarg);
			break;
		case 'c':
			contout = true;
			break;
		case 'S':
			stranded = true;
			break;
		case 'r':
			readoutput = true;
			break;
		case 'h':
			cout << helpmessage;
			return 0;
		}
	}

	string logfile_name = output_prefix + ".log";
	ofstream logfile;
	logfile.open(logfile_name);
	outlog.logfile = &logfile;

	//Parse some of the command line options to make more sense upon output to terminal and logfile.

	string minscr;
	if (score_min==0) { minscr = "UNSCORED"; }
	else { minscr = to_string(score_min); }

	string contstr;
	if (contout) { contstr = "TRUE"; }
	else { contstr = "FALSE"; }

	string crossout;
	if (crossval == 1) { crossout = "FALSE"; }
	else { crossout = to_string(crossval); }

	string strandstr;
	if (stranded) { strandstr = "TRUE"; }
	else { strandstr = "FALSE"; }

	string readoutputstr;
	if (readoutput) { readoutputstr = "TRUE"; }
	else { readoutputstr = "FALSE"; }

	outlog << "Length file: " << length_name << "\n";
	outlog << "Output prefix: " << output_prefix << "\n";
	outlog << "Iterations: " << to_string(iterations) << "\n";
	outlog << "Fixation tolerance: " << to_string(fixation) << "\n";
	outlog << "Maximum alignments: " << to_string(maxaligns) << "\n";
	outlog << "Minimum score: " << minscr << "\n";
	outlog << "Continuous output: " << contstr << "\n";
	outlog << "N-fold Cross Validation: " << crossout << "\n";
	outlog << "Strand-specific: " << strandstr << "\n";
	outlog << "Reweight fitting rate: " << to_string(fitrate) << "\n";
	outlog << "Read output: " << readoutputstr << "\n";

	//Parse length file

	ifstream length_file(length_name);

	if (!length_file) { cerr << "File " << length_name << " could not be opened for reading.\n"; return 1; }
	parseLengthFile(length_file);
	outlog << "Length file parsed.\n\n";
	length_file.close();

	//Store the list of reads files in vector rfils so it can be accessed later, multiple times if needed.

	//Create the read output file if indicated
	if (readoutput)
	{
		string readout_uniterated = output_prefix + "_reads_uniterated.bed.gz";
		read_out_unit.open(readout_uniterated.c_str());
	}

	vector< string > rfils;

	for (int fils = optind; fils < argc; fils++)
	{
		rfils.push_back(argv[fils]);
	}

	//Represents whether cross-validation is in effect
	bool cross = crossval > 1;

	//Represents whether the analysis will be "only" or "sans" the cval index
	bool onsa;

	for (int cval = 0; cval < crossval; cval++)
	{
		//Only outputs the message if crossval > 1, if there is n-fold cross-validation
		if (cross) { outlog << "Beginning cross-validation WITH ONLY section " << to_string(cval+1) << " of " << to_string(crossval) << ".\n\n"; }

		//After the first iteration, need to clear the reads and chromosome vectors as below.
		if (cval > 0)
		{
			outlog << "Clearing reads vector.\n";
			reads_vector.clear();

			outlog << "Clearing chromosome vectors.\n\n";
			for (unsigned int chromindex = 0; chromindex < tree1.size(); chromindex++)
			{
					fill(tree1[chromindex].begin(), tree1[chromindex].end(), 0);
					fill(tree2[chromindex].begin(), tree2[chromindex].end(), 0);
			}
		}

		onsa = true;
		readThroughBedgraph(rfils, iterations, fixation, output_prefix, score_min, crossval, cval, contout, cross, onsa, maxaligns);

		//No need to run through a sans analysis unless there is n-fold cross-validation
		if (cross)
		{
			outlog << "Beginning cross-validation WITHOUT section " << to_string(cval+1) << " of " << to_string(crossval) << ".\n\n";

			outlog << "Clearing reads vector.\n";
			reads_vector.clear();

			outlog << "Clearing chromosome vectors.\n\n";
			for (unsigned int chromindex = 0; chromindex < tree1.size(); chromindex++)
			{
				fill(tree1[chromindex].begin(), tree1[chromindex].end(), 0);
				fill(tree2[chromindex].begin(), tree2[chromindex].end(), 0);
			}

			onsa = false;
			readThroughBedgraph(rfils, iterations, fixation, output_prefix, score_min, crossval, cval, contout, cross, onsa, maxaligns);
		}
	}

	logfile.close();
	read_out_unit.close();
	return 0;
}

