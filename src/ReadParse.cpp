/*
 * ReadParse.cpp
 *
 *  Created on: Jun 24, 2020
 *      Author: rohan
 */

#include "BITOps.h"
#include "ReadParse.h"

using namespace std;

float probMap(int score, float score_min)
{
	if (score_min == 0) { return 1; } // Switches off probability

	float ratio = score / score_min;
	if (ratio <= 0.2) { return 0.9999369; } // MAPQ = 42
	else if (ratio <= 0.3) { return 0.9999; } // MAPQ = 40
	else if (ratio <= 0.4) { return 0.9960189; } // MAPQ = 24
	else if (ratio <= 0.5) { return 0.9949881; } // MAPQ = 23
	else if (ratio <= 0.6) { return 0.8415107; } // MAPQ = 8
	else if (ratio <= 0.7) { return 0.4988128; } // MAPQ = 3
	else { return 0; } // MAPQ = 0
}

void endRead(int& id_counter, int& not_added_fix, int numaligns, int maxaligns)
{
	int nummaps = reads_vector[id_counter].size();

	if (numaligns > maxaligns)
	{
		reads_vector.pop_back();
		id_counter--;
		return;
	}

	if (nummaps > 0)
	{
		float probsum = 0;

		for (int i = 0; i < nummaps; i++)
		{
			probsum += reads_vector[id_counter][i].prob;
		}

		if (probsum > 0)
		{
			for (int i = 0; i < nummaps; i++)
			{
				readMap tempread = reads_vector[id_counter][i];
				int chromindex = tempread.chrom;
				rangeUpdate(tree1[chromindex], tree2[chromindex], tempread.start, tempread.stop, tempread.prob / probsum, counter_to_length[chromindex]);
				reads_vector[id_counter][i].weight = tempread.prob / probsum;
			}
		}
	}

	if (nummaps == 1 || nummaps == 0)
	{
		reads_vector.pop_back();
		id_counter--;
		not_added_fix += nummaps;
	}
}

void nextRead(int& id_counter)
{
	vector<readMap> new_read_vector;
	reads_vector.push_back(new_read_vector);
	id_counter++;
}

void parseReadsFile(igzstream& reads_file, string& read_ids, int& id_counter, int& raw_counter, int crossval, int cval, int maxaligns, float score_min, bool onsa)
{
	string chrom;
	int start, stop, as, ys;
	string read_str, as_str, ys_str;
	int not_added_fix = 0;
	int not_added_prob = 0;
	int numaligns = 0;
	bool checkadd;

	while (reads_file >> chrom >> start >> stop >> read_str >> as_str >> ys_str)
	{
		if (read_str != read_ids)
		{
			checkadd = (crossval==1) || (raw_counter % crossval == cval && onsa) || (raw_counter % crossval != cval && !onsa);
			if (checkadd)
			{
				if (id_counter >= 0)
				{
					endRead(id_counter, not_added_fix, numaligns, maxaligns);
				}

				nextRead(id_counter);
				numaligns = 0;
			}
			raw_counter++;
			read_ids = read_str;
		}

		if(checkadd)
		{
			as_str.erase(0, 5);
			ys_str.erase(0, 5);
			as = stoi(as_str);
			ys = stoi(ys_str);

			float pm = probMap(as + ys, score_min);
			if (pm > 0)
			{
				readMap new_read_map = { chrom_to_counter[chrom], start, stop, 1, 0 , pm };
				reads_vector[id_counter].push_back(new_read_map);
			}
			else { not_added_prob++; }

			numaligns++;
		}
	}

	if (!checkadd) { endRead(id_counter, not_added_fix, numaligns, maxaligns); }

	outlog << "Not added due to low map quality: " << not_added_prob << "\n";
	outlog << "Not added due to fixation: " << not_added_fix << "\n";
}

void parseReadsFile(ifstream& reads_file, string& read_ids, int& id_counter, int& raw_counter, int crossval, int cval, int maxaligns, float score_min, bool onsa)
{
	string chrom;
	int start, stop, as, ys;
	string read_str, as_str, ys_str;
	int not_added_fix = 0;
	int not_added_prob = 0;
	int numaligns = 0;
	bool checkadd;

	while (reads_file >> chrom >> start >> stop >> read_str >> as_str >> ys_str)
	{
		if (read_str != read_ids)
		{
			checkadd = (crossval==1) || (raw_counter % crossval == cval && onsa) || (raw_counter % crossval != cval && !onsa);
			if (checkadd)
			{
				if (id_counter >= 0)
				{
					endRead(id_counter, not_added_fix, numaligns, maxaligns);
				}

				nextRead(id_counter);
				numaligns = 0;
			}
			raw_counter++;
			read_ids = read_str;
		}

		if(checkadd)
		{
			as_str.erase(0, 5);
			ys_str.erase(0, 5);
			as = stoi(as_str);
			ys = stoi(ys_str);

			float pm = probMap(as + ys, score_min);
			if (pm > 0)
			{
				readMap new_read_map = { chrom_to_counter[chrom], start, stop, 1, 0 , pm };
				reads_vector[id_counter].push_back(new_read_map);
			}
			else { not_added_prob++; }

			numaligns++;
		}
	}

	if (!checkadd) { endRead(id_counter, not_added_fix, numaligns, maxaligns); }

	outlog << "Not added due to low map quality: " << not_added_prob << "\n";
	outlog << "Not added due to fixation: " << not_added_fix << "\n";
}
