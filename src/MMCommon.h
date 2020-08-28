/*
 * MMCommon.h
 *
 *  Created on: Jun 24, 2020
 *      Author: rohan
 */

#ifndef MMCOMMON_H_
#define MMCOMMON_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>

using namespace std;

typedef vector< vector<float> > treesVec;
typedef vector<float> treeVec;

typedef vector< vector<double> > treesDoub;
typedef vector<double> treeDoub;

struct readMap
{
	int chrom, start, stop;
	float count, weight, prob;
	bool strand; // TRUE if plus, FALSE if minus
};

struct outStream
{
	ofstream* logfile;

	template<typename T> outStream& operator<<(const T& val)
	{
		cout << val;
		*logfile << val;
		return *this;
	}
};


//Declaration of global variables, initialized in Multimap.cpp

extern outStream outlog;
extern treesVec tree1, tree1neg;
extern treesDoub tree2, tree2neg;
extern vector< vector<readMap> > reads_vector;
extern map<string, int> chrom_to_counter;
extern vector<string> counter_to_chrom;
extern vector<int> counter_to_length;
extern bool stranded;

#endif /* MMCOMMON_H_ */
