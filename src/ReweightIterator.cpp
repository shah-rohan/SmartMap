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

void removeFixedRead(int i)
{
	if (i != int(reads_vector.size()) - 1)
	{
		reads_vector[i] = std::move(reads_vector.back());
	}
	reads_vector.pop_back();
}

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
			bool fixed = true;

			if (reads_vector[i].size() == 1)
			{
				removeFixedRead(i);
				fixremove++;
				continue;
			}

			int nummaps = reads_vector[i].size();
			float countsum = 0;

			for (int j = 0; j < nummaps; j++)
			{
				int chromindex = reads_vector[i][j].chrom;
				reads_vector[i][j].count = rangeSum(tree1[chromindex], tree2[chromindex], reads_vector[i][j].start, reads_vector[i][j].stop);
			}

			for (int j = 0; j < nummaps; j++)
			{
				if (reads_vector[i][j].count * reads_vector[i][j].prob < 0)
				{
					continue;
				}
				countsum += reads_vector[i][j].count * reads_vector[i][j].prob * 100.0 / (reads_vector[i][j].stop - reads_vector[i][j].start);
			}

			for (int j = 0; j < nummaps; j++)
			{
				readMap tempread = reads_vector[i][j];
				float new_weight = (tempread.count * tempread.prob * 100.0 / (tempread.stop - tempread.start)) / countsum;
				if (tempread.count * tempread.prob < 0 || countsum <= 0) { new_weight = 0; }
				float dw = new_weight - tempread.weight;

				if (dw == 0)
				{
					continue;
				}

				int chromindex = tempread.chrom;
				rangeUpdate(tree1[chromindex], tree2[chromindex], tempread.start, tempread.stop, dw, counter_to_length[chromindex]);
				reads_vector[i][j].weight = new_weight;
				if (new_weight < (1-fixation) * tempread.weight || new_weight > (1+fixation) * tempread.weight)
				{
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
