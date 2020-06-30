/*
 * BITOps.cpp
 *
 *  Created on: Jun 24, 2020
 *      Author: rohan
 */

#include "BITOps.h"

using namespace std;

#define lsb(i) ((i) & -(i))

float treeSum(const treeVec& tree, int index)
{
	float sum = 0;
	while (index > 0)
	{
		sum += tree[index];
		index -= lsb(index);
	}

	return sum;
}

double treeSum(const treeDoub& tree, int index)
{
	float sum = 0;
	while (index > 0)
	{
		sum += tree[index];
		index -= lsb(index);
	}

	return sum;
}

void treeUpdate(treeVec& tree, int index, float val, int tree_size)
{
	while (index <= tree_size)
	{
		tree[index] += val;
		index += lsb(index);
	}
}

void treeUpdate(treeDoub& tree, int index, float val, int tree_size)
{
	while (index <= tree_size)
	{
		tree[index] += val;
		index += lsb(index);
	}
}

double pointSum(const treeVec& tree1, const treeDoub& tree2, int index)
{
	return double(treeSum(tree1, index)) * index - treeSum(tree2, index);
}

//Returns the range from the interval [start, stop) = [start, stop - 1]
float rangeSum(const treeVec& tree1, const treeDoub& tree2, int start, int stop)
{
	return float(pointSum(tree1, tree2, stop - 1) - pointSum(tree1, tree2, start - 1));
}

//Updates the range from the interval [start, stop) = [start, stop - 1]
void rangeUpdate(treeVec& tree1, treeDoub& tree2, int start, int stop, float val, int tree_size)
{
	treeUpdate(tree1, start, val, tree_size);
	treeUpdate(tree1, stop, -val, tree_size);
	treeUpdate(tree2, start, val * (start - 1), tree_size);
	treeUpdate(tree2, stop, -val * (stop - 1), tree_size);
}


