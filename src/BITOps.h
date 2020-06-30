/*
 * BITOps.h
 *
 *  Created on: Jun 24, 2020
 *      Author: rohan
 */

#ifndef BITOPS_H_
#define BITOPS_H_

#include <vector>
#include "MMCommon.h"

using namespace std;

float treeSum(const treeVec& tree, int index);
float rangeSum(const treeVec& tree1, const treeDoub& tree2, int start, int stop);
void rangeUpdate(treeVec& tree1, treeDoub& tree2, int start, int stop, float val, int tree_size);

#endif /* BITOPS_H_ */
