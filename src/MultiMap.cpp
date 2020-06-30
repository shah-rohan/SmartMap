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

#include "gzstream.h"
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
treesVec tree1;
treesDoub tree2;
vector< vector<readMap> > reads_vector;
map<string, int> chrom_to_counter;
vector<string> counter_to_chrom;
vector<int> counter_to_length;

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

void readThroughBedgraph(const vector<string> rfils, int iterations, float fixation, string output_prefix, float score_min, int crossval, int cval, bool contout, bool cross, bool onsa, int maxaligns)
{
	string read_ids = "";
	int id_counter = -1;
	int raw_counter = 0;

	for (unsigned int fils = 0; fils < rfils.size(); fils++)
	{
		outlog << "Now parsing read file " << rfils[fils] << "\n";
		auto start = high_resolution_clock::now();

		bool checkgz = hasEnding(rfils[fils], ".gz");

		if (checkgz)
		{
			igzstream reads_file(rfils[fils].c_str());
			if (!reads_file) { cerr << "File " << rfils[fils] << " could not be opened for reading."; return; }
			parseReadsFile(reads_file, read_ids, id_counter, raw_counter, crossval, cval, maxaligns, score_min, onsa);
			reads_file.close();
		}
		else
		{
			ifstream reads_file(rfils[fils]);
			if (!reads_file) { cerr << "File " << rfils[fils] << " could not be opened for reading."; return; }
			parseReadsFile(reads_file, read_ids, id_counter, raw_counter, crossval, cval, maxaligns, score_min, onsa);
			reads_file.close();
		}

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<seconds>(stop - start);
		outlog << "Finished parsing read file " << rfils[fils] << " in " << duration.count() << " seconds\n\n";
	}
	read_ids.clear();

	string reweight_prefix = output_prefix;
	if (cross)
	{
		if (onsa) { reweight_prefix = reweight_prefix + "_" + to_string(crossval) + "-fold-only-" + to_string(cval+1); }
		else { reweight_prefix = reweight_prefix + "_" + to_string(crossval) + "-fold-sans-" + to_string(cval+1); }
	}

	reweightIterator(iterations, fixation, reweight_prefix, contout);

	bool finalwrite = iterations < 10 || iterations % 10 == 0;

	if (!contout || !finalwrite)
	{
		string filename = reweight_prefix + ".bedgraph.gz";
		outlog << "Writing bedgraph to " << filename << "\n";
		writeBedgraphOutput(filename);
	}
}

int main(int argc, char* argv[])
{
	int opt;
	int iterations = 1, maxaligns = 50, crossval = 1;
	float fixation = 0, score_min = 0;
	string length_name, output_prefix;
	bool contout = false;

	string helpmessage = "MultiMap for analysis of ambiguously mapping reads in ChIP-seq.\n\n"
			"Usage: MultiMap [options] [bed or bed.gz file input(s)]\n\n"
			"Required options:\n"
			"-g : Genome length file listing all chromosomes and lengths. Chromosomes will appear in this order in output file.\n"
			"-o : Output prefix used for output bedgraph and log files.\n\n"
			"Non-required options:\n"
			"-i : Number of iterations. Default 1.\n"
			"-x : Fixation percentage threshold. Default 0.\n"
			"-m : Maximum number of alignments for a read to be processed. Default 50.\n"
			"-s : Minimum score for Bowtie2 display. Default 0 (unscored).\n"
			"-v : N for N-fold cross-validation. Default 1 (no cross-validation).\n"
			"-c : Flag for continuous output bedgraphs. Default off.\n"
			"-h : Display help message.\n\n"
			"Developed by Rohan Shah (rohanshah@uchicago.edu).\n";

	//Parse command options

	while ((opt = getopt(argc, argv, "g:o:i:x:m:s:v:ch")) != -1)
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
		case 'c':
			contout = true;
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

	string minscr;
	if (score_min==0) { minscr = "UNSCORED"; }
	else { minscr = to_string(score_min); }

	string contstr;
	if (contout) { contstr = "TRUE"; }
	else { contstr = "FALSE"; }

	string crossout;
	if (crossval == 1) { crossout = "FALSE"; }
	else { crossout = to_string(crossval); }

	outlog << "Length file: " << length_name << "\n";
	outlog << "Output prefix: " << output_prefix << "\n";
	outlog << "Iterations: " << to_string(iterations) << "\n";
	outlog << "Fixation tolerance: " << to_string(fixation) << "\n";
	outlog << "Maximum alignments: " << to_string(maxaligns) << "\n";
	outlog << "Minimum score: " << minscr << "\n";
	outlog << "Continuous output: " << contstr << "\n";
	outlog << "N-fold Cross Validation: " << crossout << "\n";

	//Parse length file

	ifstream length_file(length_name);

	if (!length_file) { cerr << "File " << length_name << " could not be opened for reading.\n"; return 1; }
	parseLengthFile(length_file);
	outlog << "Length file parsed.\n\n";
	length_file.close();

	//Parse the reads file

	vector< string > rfils;

	for (int fils = optind; fils < argc; fils++)
	{
		rfils.push_back(argv[fils]);
	}

	bool cross = crossval > 1;
	bool onsa;

	for (int cval = 0; cval < crossval; cval++)
	{
		if (cross) { outlog << "Beginning cross-validation WITH ONLY section " << to_string(cval+1) << " of " << to_string(crossval) << ".\n\n"; }

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
	return 0;
}

