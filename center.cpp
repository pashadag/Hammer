#include<cassert>
#include<cmath>
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<cstdarg>
#include<algorithm>

using namespace std;
ostringstream oMsg;
string sbuf;
#include "defs.h"


double errorRate = 0.01;
double ambigConst = 1;
int nthreads;
string ufFilename;
double threshold;
int savelimit;

class Read {
public:
	string id;
	string seq;
	int count;
	float freq;
	bool operator<(Read x) const {
		return x.seq < seq;
	}
	bool operator==(string x) const {
		return x == this->seq;
	}
};

/*bool findPseudoKmers (string x, string y, string & new1, string & new2) {
	assert(x.length() == y.length());
	int diff = 0;
	char x1, x2, y1, y2;
	int pos1, pos2;
	for (int i = 0; i < x.length(); i++) {
		if (x[i] == y[i]) continue;
		diff++;
		if (diff == 1) {
			pos1 = i;
			x1 = x[i];
			y1 = y[i];
		} else if (diff == 2) {
			pos2 = i;
			x2 = x[i];
			y2 = y[i];
		} else {
			return false;
		}
	}
	if (diff == 2) {
		new1 = x;
		new1[pos1] = y1;
		new2 = y;
		new2[pos2] = x2;
		return true;
	}
	return false;
}
*/

double calcMultCoef(vector<int> & distances, vector<Read> & kmers) {
	int kmersize = kmers[0].seq.size();
	double prob = 0;
	double theta;
	//cout << "calcMultCoef.  distances = " << distances << "\t\t\t"; for (int i = 0; i < kmers.size(); i++) cout << kmers[i].mult << "\t";
	//cout << "\nadding:\t";
	for (size_t i = 0; i < kmers.size(); i++) {
		theta = kmers[i].count * -((kmersize - distances[i]) * log(1 - errorRate) + distances[i] * log(errorRate));
		//cout << theta << "\t";
		prob = prob + theta;
	}
	//cout << endl;

	return prob;
}

string find_consensus(vector<Read> & block) {
	string c;
	
	for (int i = 0; i < block[0].seq.length(); i++) {
		int scores[4] = {0,0,0,0};
		for (int j = 0; j < block.size(); j++) {
			scores[nt2num(block[j].seq.at(i))]++;
		}
		c.push_back(num2nt(argmax(scores, 4)));
	}
	return c;
}

		/*
		for (int i = 0; i < origBlockSize; i++) {
			for (int j = 0; j < origBlockSize; j++) {
				string new1, new2;
				if (findPseudoKmers(block[i].seq, block[j].seq, new1, new2)) {
					Read cur;
					new1 = rcnorm(new1);
					new2 = rcnorm(new2);
					vector<Read>::iterator it1 = find(block.begin(), block.end(), new1);
					vector<Read>::iterator it2 = find(block.begin(), block.end(), new2);
					if (it1 == block.end()) {
						cur.seq = new1;
						cur.count = 0;
						block.push_back(cur);
					} 
					if (it2 == block.end()) {
						cur.seq = new2;
						cur.count = 0;
						block.push_back(cur);
					} 
				}
			}
		}
		*/

void process_block(vector<Read> & block, string blockNum, double threshold, ofstream & outf) {
	int origBlockSize = block.size();
	if (block.size() == 0) return;
	vector<double> multiCoef(block.size() + 1,1000000);//add one for the consensus
	vector<int> distance(block.size() + 1 , 0);  //add one for the consensus
	vector< vector<int> > distances(block.size() + 1 , distance); //add one for the consensus
	//sort (block.begin(), block.end());
	string bestkmer;
	bool change = false;
	string reason = "noreason";
	bool bigboy = false;

	if (block.size() > 1) {

		//This is a shortcut that deviates from the definition of the algorithm but speeds everything up
		//if one node has a count more than double of another one
		int max1 = 0;
		int max2 = 1;
		if (block[max1].count < block[max2].count) swap(max1, max2);
		for (int i = 2; i < block.size(); i++) {
			if (block[i].count > block[max1].count) {
				max2 = max1;
				max1 = i;
			} else if (block[i].count > block[max2].count) {
				max2 = i;
			}
		}

		if (block[max1].count >= 2* block[max2].count && block[max1].count >= 20) {
			change = true;
			bestkmer = block[max1].seq;
			bigboy = true;
		} else {
			//Add a consensus sequence
			Read consensus;
			consensus.seq = find_consensus(block);
			consensus.count = 0;
			bool found = false;
			for (int i = 0; i < block.size(); i++) {
				if (block[i].seq == consensus.seq) {
					found = true;
					break;
				}
			}
			if (!found) block.push_back(consensus);
			

			//Calculate distance matrix
			for (int i = 0; i < block.size(); i++) {
				distances[i][i] = 0;
				for (int j = i + 1; j < block.size(); j++) {
					distances[i][j] = hamdist(block[i].seq, block[j].seq, ANY_STRAND);
					distances[j][i] = distances[i][j];
				}
			}

			//Calculate multinmial coefficient
			for (int i = 0; i < block.size(); i++) {
				multiCoef[i] = calcMultCoef(distances[i], block);
			}

			//Pick a center node
			double low1 = 1000000000;
			int ind1 = -1;
			int ind2 = -1;
			double low2 = 1000000000;
			for (int i = 0; i < block.size(); i++) {
				if (multiCoef[i] < low1) {
					low2 = low1;
					ind2 = ind1;
					low1 = multiCoef[i];
					ind1 = i;
				} else if (multiCoef[i] < low2) {
					low2 = multiCoef[i];
					ind2 = i;
				}
			}
			if (low1 < low2 - ambigConst) {
				change = true; 
				if (ind1 == origBlockSize)  reason = "pseudo";  //consensus was picked
			}
			bestkmer = block[ind1].seq;
		}
	} 

	for (int i = 0; i < origBlockSize; i++) {
		string replacement = bestkmer;
		if (block.size() == 1) { //singleton
			if (block[i].freq >= threshold) { //keep reads
				replacement = block[i].seq;
				reason =  "goodSingleton";
				change = false;
			} else {
				replacement = "removed";
				reason =  "badSingleton";
				change = true;
			}
		} else { //multiton
			if (change) {
				if (reason != "pseudo") {
					if (bestkmer == block[i].seq) {
						reason = "center";
					} else if (block[i].freq >= savelimit) { //save from destruction cause its too frequent
						replacement = block[i].seq;
						reason = "saved";
					} else {
						reason = "change";
					} 
				}
			} else { //ambig, treat like singleton
				if (block[i].freq >= threshold) { //keep reads
					replacement =  block[i].seq;
					reason =  "goodAmbig";
				} else {
					replacement = "removed";
					reason =  "badAmbig";
				}
			}
		}

		outf << blockNum << "\t" << block[i].seq << "\t" << replacement << "\t";
		outf << block[i].count << "\t" << block[i].freq << "\t" << block.size();
		outf <<  "\t" << multiCoef[i] << "\t" << reason << "\t";
		outf << distances[i][0]; 
		if( !bigboy) {  
			for (int j = 1; j < distances[i].size(); j++) outf << "_" << distances[i][j];
		}
		outf << endl;
	}


	outf << endl;
	return;
}

void * onethread(void * params) {
	int thread = *((int *) params);

	ifstream inf;
	open_file(inf, ufFilename);
	ofstream outf;
	open_file(outf, "reads.uf.corr." + make_string(thread));

	vector<Read> block;
	int counter=0;
	vector<string> row;
	string curBlockNum;
	string lastBlockNum = "GO LAKERS!";
	while (get_row_whitespace(inf, row)) {
		if (row.size() < 1 || row[0] != "ITEM") {
			//outf << row << endl;
			continue;
		}
		if (atoi(row[1].c_str()) % nthreads != thread) continue;
		if (++counter % 1000000 == 0) cerr << "Processed " << add_commas(counter) << ", ";
		Read cur;
		curBlockNum = row[1];
		cur.id = row[2];
		cur.seq = row[3];
		cur.count = atoi(row[4].c_str());
		cur.freq = atof(row[5].c_str());
		if (lastBlockNum == curBlockNum) { //add to current reads
			block.push_back(cur);
		} else {
			process_block(block,lastBlockNum, threshold, outf);
			block.clear();
			block.push_back(cur);
		}
		lastBlockNum = curBlockNum;
	}
	process_block(block,curBlockNum, threshold, outf);
	cerr << "Finished\n";
	inf.close();
	outf.close();
	pthread_exit(NULL);


}
int main(int argc, char * argv[]) {
	assert(argc == 5);

	ufFilename = argv[1];
	threshold = atof(argv[2]);
	nthreads = atoi(argv[3]);
	savelimit = atoi(argv[4]);
	/*
	   if (threshold == -1) {
	   cerr << "First phase...\n";
	   threshold = get_threshold(ufFilename);
	   }
	 */
	cerr << "Second phase (threshold = " << threshold << ")...\n";


	//start threads
	pthread_t thread[nthreads];
	int params[nthreads];
	for (int i = 0; i < nthreads; i++) {
		params[i] = i;
		pthread_create(&thread[i], NULL, onethread, (void *) &params[i]);
	}
	for (int i = 0; i < nthreads; i++) {
		pthread_join(thread[i], NULL);
	}

	return 0;
}

/*

   double get_min(vector<double> & counts) {
//cout << endl << counts << "\tCOUNTS" << endl;
//processes the counts
int res = 10;
vector<double> bigbins(100,0);
vector<double> smallbins(100*res,0);

//put them into big and small bins first.
for (int i = 0; i < counts.size(); i++) {
int bin = int(counts[i]);
if (bin < 100) bigbins.at(bin) += 1;
bin = int(counts[i] * res);
if (bin < 100 * res ) smallbins.at(bin) += 1;
}


ofstream dump;
open_file(dump, "singletonCounts.txt");
dump << bigbins << "\tBIG_BINS" << endl; 
dump << smallbins << "\tSMALL_BINS" << endl; 
dump.close();

//find first min big bin
int minBigBin = -1;
for (int i = 1; i < bigbins.size() - 1; i++) {
if (bigbins[i-1] > bigbins[i] && bigbins[i] < bigbins[i+1]) {
minBigBin = i;
break;
}
}
//cout << bigbins << "\tBIG_BINS" << endl; cout << minBigBin << "\tminBigBin" << endl;

if (minBigBin == -1) {
cerr << "Cannot find minimum in frequency distribution!  Will use 0.\n";
return 0;
}


//find first min bin
int minSmallBin = -1;
for (int i = max(0, (minBigBin -1 ) * res) + 1; i < min(int(smallbins.size()), res * (minBigBin + 1) ) - 1; i++) {
if (smallbins[i-1] > smallbins[i] && smallbins[i] < smallbins[i+1]) {
minSmallBin = i;
break;
}
}


return minBigBin;

double retval = minBigBin - 1 + (minSmallBin / double(res));
assert (minSmallBin >= 0);
return retval;
}

double get_total_counts(vector<Read> & block) {

// returns the sum of the weights in this block, -1 if its not a singleton
//Create list of kmers and their multiplicities
if (block.size() == 0) return -1;

//vector<Read> orBlock(block);
double sum = block[0].count;
string seq = rcnorm(block[0].seq);
for (int i = 1; i < block.size(); i++) {
if (rcnorm(block[i].seq) != seq) return -1;
sum += block[i].count;
}
return sum;
}


double get_threshold(string ufFilename) {
	double retval;
	ifstream inf;
	open_file(inf, ufFilename);


	vector<double> counts;
	vector<Read> block;
	int counter=0;
	double totCount;
	vector<string> row;
	string curBlockNum;
	string lastBlockNum = "GO LAKERS!";
	while (get_row_whitespace(inf, row)) {
		if (++counter % 1000000 == 0) cerr << "Processed " << add_commas(counter) << ", ";
		if (row.size() < 1 || row[0] != "ITEM") {
			continue;
		}
		Read cur;
		curBlockNum = row[1];
		cur.id = row[2];
		cur.seq = row[3];
		cur.count = atof(row[4].c_str());

		if (lastBlockNum == curBlockNum) { //add to current reads
			block.push_back(cur);
		} else {
			totCount = get_total_counts(block);
			if (totCount >= 0) counts.push_back(totCount);
			block.clear();
			block.push_back(cur);
		}
		lastBlockNum = curBlockNum;
	}
	cerr << endl;
	totCount = get_total_counts(block);
	if (totCount >= 0) counts.push_back(totCount);
	retval = get_min(counts);
	inf.close();
	return retval;

}

*/
