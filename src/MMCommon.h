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
extern treesVec tree1;
extern treesDoub tree2;
extern vector< vector<readMap> > reads_vector;
extern map<string, int> chrom_to_counter;
extern vector<string> counter_to_chrom;
extern vector<int> counter_to_length;

#endif /* MMCOMMON_H_ */
