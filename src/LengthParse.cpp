/*
 * LengthParse.cpp
 *
 *  Created on: Jun 24, 2020
 *      Author: rohan
 */

#include "LengthParse.h"

using namespace std;

//Parse the length file into map and vector to associate chromosome to counter, chromosome lengths, and add to the trees
void parseLengthFile(ifstream& length_file)
{
	string chrom;
	int length;
	int counter = 0;

	if (!length_file)
	{
		cerr << "Can't open " << endl;
	}

	while (length_file >> chrom >> length)
	{
		outlog << "Now parsing chromosome: " << chrom << "\n";
		chrom_to_counter[chrom] = counter;
		counter_to_chrom.push_back(chrom);
		counter_to_length.push_back(length);
		treeVec chromtemp(length + 1, 0);
		treeDoub chromtemp2(length + 1, 0);
		tree1.push_back(chromtemp);
		tree2.push_back(chromtemp2);
		counter++;
	}
}
