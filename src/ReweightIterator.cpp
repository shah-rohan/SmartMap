/*
 * ReweightIterator.cpp
 *
 *  Created on: Jun 24, 2020
 *      Author: rohan
 */
#include <chrono>

#include "BedgraphWrite.h"
#include "BITOps.h"
#include "ReweightIterator.h"

using namespace std;
using namespace std::chrono;

//Removes a read from reads_vector with index i; move operator to be avoided if it's the last read already to avoid move to self.
void removeFixedRead(int i)
{
	if (i != int(reads_vector.size()) - 1)
	{
		reads_vector[i] = std::move(reads_vector.back());
	}
	reads_vector.pop_back();
}

//Runs through the reads_vector and iteratively reassigns weights as per the posterior probability algorithm.
void reweightIterator(int iterations, float fixation, string output_prefix, bool contout)
{
	for (int itercounter = 1; itercounter <= iterations; itercounter++)
	{
		outlog << "Begin iteration " << itercounter << "\n";
		int fixremove = 0;
		auto start = high_resolution_clock::now();
		outlog << "Total reads: " << to_string(reads_vector.size()) << "\n";

		for (int i = reads_vector.size()-1; i >= 0; i--)
		{
			//Assumes at the beginning of each read that it has reached fixation until proven otherwise
			bool fixed = true;

			//Really should be unnecessary due to the removal of fixed reads in ReadParse, but nothing wrong with an extra check.
			if (reads_vector[i].size() == 1)
			{
				removeFixedRead(i);
				fixremove++;
				continue;
			}

			int nummaps = reads_vector[i].size();
			float countsum = 0;

			//Updates the counts of each map from the trees.
			for (int j = 0; j < nummaps; j++)
			{
				int chromindex = reads_vector[i][j].chrom;
				reads_vector[i][j].count = rangeSum(tree1[chromindex], tree2[chromindex], reads_vector[i][j].start, reads_vector[i][j].stop);
			}

			//Adds up the counts for the reads
			for (int j = 0; j < nummaps; j++)
			{
				if (reads_vector[i][j].count * reads_vector[i][j].prob < 0)
				{
					continue;
				}
				countsum += reads_vector[i][j].count * reads_vector[i][j].prob * 100.0 / (reads_vector[i][j].stop - reads_vector[i][j].start);
			}

			//Reassign the weights based on the probability of proper alignment, counts, and the countsum
			for (int j = 0; j < nummaps; j++)
			{
				readMap tempread = reads_vector[i][j];
				float new_weight = (tempread.count * tempread.prob * 100.0 / (tempread.stop - tempread.start)) / countsum;
				if (tempread.count * tempread.prob < 0 || countsum <= 0) { new_weight = 0; }
				float dw = new_weight - tempread.weight;

				//If there's no changes in the map weight, then there's no need to go through the motions of adding zero to the trees.
				if (dw == 0)
				{
					continue;
				}

				int chromindex = tempread.chrom;
				rangeUpdate(tree1[chromindex], tree2[chromindex], tempread.start, tempread.stop, dw, counter_to_length[chromindex]);
				reads_vector[i][j].weight = new_weight;
				if (new_weight < (1-fixation) * tempread.weight || new_weight > (1+fixation) * tempread.weight)
				{
					//This part is only reached if there is a change dw outside the fixation tolerance (by default: 0). One unfixed map unfixes the read.
					fixed = false;
				}
			}

			if (fixed)
			{
				removeFixedRead(i);
				fixremove++;
			}
		}

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<seconds>(stop - start);

		outlog << "Reads removed due to fixation: " << fixremove << "\n";
		outlog << "Completed iteration " << itercounter << " in " << duration.count() << " seconds\n\n";

		//Continuous output bedgraph option: prints for the first 10 iterations and every 10th iteration thereafter
		if (contout)
		{
			if (itercounter < 10 || itercounter % 10 == 0)
			{
				outlog << "Writing bedgraph for iteration " << itercounter << "\n\n";
				string filename = output_prefix + "_iteration-" + to_string(itercounter) + ".bedgraph.gz";
				writeBedgraphOutput(filename);
			}
		}
	}
}
