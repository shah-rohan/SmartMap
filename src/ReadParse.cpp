/*
 * ReadParse.cpp
 *
 *  Created on: Jun 24, 2020
 *      Author: rohan
 */

#include "BITOps.h"
#include "ReadParse.h"

using namespace std;

//Assigns probability of correct alignment as per the AS: and YS: scores from Bowtie2 as if the map was a uniread based on MAPQ.
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

//Software for the end of a read
void endRead(int& id_counter, int& not_added_fix, int numaligns, int maxaligns)
{
	int nummaps = reads_vector[id_counter].size();

	//If the number of maps processed in that read is greater than the number of maximum allowable number of maps, remove the read and decrement id_counter
	if (numaligns > maxaligns)
	{
		reads_vector.pop_back();
		id_counter--;
		return;
	}

	if (nummaps > 0)
	{
		float probsum = 0;

		//Sum of the probabilities of the maps for that read
		for (int i = 0; i < nummaps; i++)
		{
			probsum += reads_vector[id_counter][i].prob;
		}

		if (probsum > 0)
		{
			//For each of those maps, assign the weight as the probability divided by the sum of probabilities over the read and  add to the trees.
			for (int i = 0; i < nummaps; i++)
			{
				readMap tempread = reads_vector[id_counter][i];
				int chromindex = tempread.chrom;
				if (tempread.strand)
				{
					rangeUpdate(tree1[chromindex], tree2[chromindex], tempread.start, tempread.stop, tempread.prob / probsum, counter_to_length[chromindex]);
				}
				else
				{
					rangeUpdate(tree1neg[chromindex], tree2neg[chromindex], tempread.start, tempread.stop, tempread.prob / probsum, counter_to_length[chromindex]);
				}
				reads_vector[id_counter][i].weight = tempread.prob / probsum;
			}
		}
	}

	//Fixed reads: Don't need them in the array, so just remove the read and decrement id_counter. Add 1 to fixed number if there was 1 map, 0 if 0 maps.
	if (nummaps == 1 || nummaps == 0)
	{
		reads_vector.pop_back();
		id_counter--;
		not_added_fix += nummaps;
	}
}

//Add a new vector to the reads_vector and increment id_counter to analyze the next read.
void nextRead(int& id_counter)
{
	vector<readMap> new_read_vector;
	reads_vector.push_back(new_read_vector);
	id_counter++;
}

//Read parser function for compressed files.
void parseReadsFile(igzstream& reads_file, string& read_ids, int& id_counter, int& raw_counter, int crossval, int cval, int maxaligns, float score_min, bool onsa)
{
	string chrom;
	int start, stop, as, ys;
	string read_str, as_str, ys_str;
	int not_added_fix = 0;
	int not_added_prob = 0;
	int numaligns = 0;
	/* Cross-validation functionality is sufficiently complex that it deserves a dedicated explanation here.
	 * To determine whether a read needs to be added, checkadd is computed at the detection of a read ID for the beginning of the next read (or first read).
	 * checkadd is true under one of the three conditions: no cross-validation, modulus matching cval in "only" mode, or modulus not matching cval in "sans" mode.
	 * If checkadd was true for the previous read, then the previous read is added to the chromosome trees and a new read is added.
	 * If checkadd is not true, then there's no need to add a new read vector, and the previous read will be processed with the next read for which checkadd is true.
	 * Either way, the raw_counter is incremented. For this read, values are only added to the latest vector if checkadd is true.
	 * If checkadd was false at the end of the while loop (end of reads), then the last read didn't get added to the trees, so it is added with endRead after loop. */

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
			//Remove the AS:i: and the YS:i: from the scores, respectively.
			as_str.erase(0, 5);
			ys_str.erase(0, 5);
			as = stoi(as_str);
			ys = stoi(ys_str);

			//Check the probability of the map being accurately aligned; if zero, then it will always have weight 0, so just remove the map.
			float pm = probMap(as + ys, score_min);
			if (pm > 0)
			{
				readMap new_read_map = { chrom_to_counter[chrom], start, stop, 1, 0, pm, true };
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

//Read parser for uncompressed files
void parseReadsFile(ifstream& reads_file, string& read_ids, int& id_counter, int& raw_counter, int crossval, int cval, int maxaligns, float score_min, bool onsa)
{
	string chrom;
	int start, stop, as, ys;
	string read_str, as_str, ys_str;
	int not_added_fix = 0;
	int not_added_prob = 0;
	int numaligns = 0;
		/* Cross-validation functionality is sufficiently complex that it deserves a dedicated explanation here.
		 * To determine whether a read needs to be added, checkadd is computed at the detection of a read ID for the beginning of the next read (or first read).
		 * checkadd is true under one of the three conditions: no cross-validation, modulus matching cval in "only" mode, or modulus not matching cval in "sans" mode.
		 * If checkadd was true for the previous read, then the previous read is added to the chromosome trees and a new read is added.
		 * If checkadd is not true, then there's no need to add a new read vector, and the previous read will be processed with the next read for which checkadd is true.
		 * Either way, the raw_counter is incremented. For this read, values are only added to the latest vector if checkadd is true.
		 * If checkadd was false at the end of the while loop (end of reads), then the last read didn't get added to the trees, so it is added with endRead after loop. */

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
			//Remove the AS:i: and the YS:i: from the scores, respectively.
			as_str.erase(0, 5);
			ys_str.erase(0, 5);
			as = stoi(as_str);
			ys = stoi(ys_str);

			//Check the probability of the map being accurately aligned; if zero, then it will always have weight 0, so just remove the map.
			float pm = probMap(as + ys, score_min);
			if (pm > 0)
			{
				readMap new_read_map = { chrom_to_counter[chrom], start, stop, 1, 0, pm , true};
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

void parseReadsFileStranded(igzstream& reads_file, string& read_ids, int& id_counter, int& raw_counter, int crossval, int cval, int maxaligns, float score_min, bool onsa)
{
	string chrom;
	int start, stop, as, ys;
	string read_str, as_str, ys_str, strand_str;
	int not_added_fix = 0;
	int not_added_prob = 0;
	int numaligns = 0;
		/* Cross-validation functionality is sufficiently complex that it deserves a dedicated explanation here.
		 * To determine whether a read needs to be added, checkadd is computed at the detection of a read ID for the beginning of the next read (or first read).
		 * checkadd is true under one of the three conditions: no cross-validation, modulus matching cval in "only" mode, or modulus not matching cval in "sans" mode.
		 * If checkadd was true for the previous read, then the previous read is added to the chromosome trees and a new read is added.
		 * If checkadd is not true, then there's no need to add a new read vector, and the previous read will be processed with the next read for which checkadd is true.
		 * Either way, the raw_counter is incremented. For this read, values are only added to the latest vector if checkadd is true.
		 * If checkadd was false at the end of the while loop (end of reads), then the last read didn't get added to the trees, so it is added with endRead after loop. */

	bool checkadd;

	while (reads_file >> chrom >> start >> stop >> read_str >> strand_str >> as_str >> ys_str)
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
			//Remove the AS:i: and the YS:i: from the scores, respectively.
			as_str.erase(0, 5);
			ys_str.erase(0, 5);
			as = stoi(as_str);
			ys = stoi(ys_str);
			bool strand;

			if (strand_str == "+") { strand = true; }
			else { strand = false; }

			//Check the probability of the map being accurately aligned; if zero, then it will always have weight 0, so just remove the map.
			float pm = probMap(as + ys, score_min);
			if (pm > 0)
			{
				readMap new_read_map = { chrom_to_counter[chrom], start, stop, 1, 0, pm , strand};
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

void parseReadsFileStranded(ifstream& reads_file, string& read_ids, int& id_counter, int& raw_counter, int crossval, int cval, int maxaligns, float score_min, bool onsa)
{
	string chrom;
	int start, stop, as, ys;
	string read_str, as_str, ys_str, strand_str;
	int not_added_fix = 0;
	int not_added_prob = 0;
	int numaligns = 0;
		/* Cross-validation functionality is sufficiently complex that it deserves a dedicated explanation here.
		 * To determine whether a read needs to be added, checkadd is computed at the detection of a read ID for the beginning of the next read (or first read).
		 * checkadd is true under one of the three conditions: no cross-validation, modulus matching cval in "only" mode, or modulus not matching cval in "sans" mode.
		 * If checkadd was true for the previous read, then the previous read is added to the chromosome trees and a new read is added.
		 * If checkadd is not true, then there's no need to add a new read vector, and the previous read will be processed with the next read for which checkadd is true.
		 * Either way, the raw_counter is incremented. For this read, values are only added to the latest vector if checkadd is true.
		 * If checkadd was false at the end of the while loop (end of reads), then the last read didn't get added to the trees, so it is added with endRead after loop. */

	bool checkadd;

	while (reads_file >> chrom >> start >> stop >> read_str >> strand_str >> as_str >> ys_str)
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
			//Remove the AS:i: and the YS:i: from the scores, respectively.
			as_str.erase(0, 5);
			ys_str.erase(0, 5);
			as = stoi(as_str);
			ys = stoi(ys_str);
			bool strand;

			if (strand_str == "+") { strand = true; }
			else { strand = false; }

			//Check the probability of the map being accurately aligned; if zero, then it will always have weight 0, so just remove the map.
			float pm = probMap(as + ys, score_min);
			if (pm > 0)
			{
				readMap new_read_map = { chrom_to_counter[chrom], start, stop, 1, 0, pm , strand};
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
