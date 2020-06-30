/*
 * ReadParse.h
 *
 *  Created on: Jun 24, 2020
 *      Author: rohan
 */

#ifndef READPARSE_H_
#define READPARSE_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <iomanip>
#include <chrono>

#include "gzstream.h"
#include "MMCommon.h"

using namespace std;
using namespace std::chrono;

void parseReadsFile(igzstream& reads_file, string& read_ids, int& id_counter, int& raw_counter, int crossval, int cval, int maxaligns, float score_min, bool onsa);
void parseReadsFile(ifstream& reads_file, string& read_ids, int& id_counter, int& raw_counter, int crossval, int cval, int maxaligns, float score_min, bool onsa);

#endif /* READPARSE_H_ */
