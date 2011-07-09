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
#include<cstdarg>
#include<algorithm>

using namespace std;
#include "defs.h"
#include "sa.h"

struct compeqNum2NumFunctor {

	bool operator() (long x, long y) const {
		string xstr = genome.substr(x,kmersize);
		string ystr = genome.substr(y,kmersize);
		return rcnorm(xstr) == rcnorm(ystr);
	}

	string genome;
	int kmersize;
};

struct compltNum2NumFunctor {
	bool operator() (long x, long y) const {
		if (x==y) return false;
		bool flipx = false;  // are we revcomping x?
		bool flipy = false;
		bool detx = false; //have we determine whether we are revcomping x?
		bool dety = false;
		for (int i = 0; i < kmersize; i++) {
			char xchar, ychar;
			if (detx && !flipx) {
				xchar = genome.at(x+i);
			} else if (detx && flipx) {
				xchar = revcomp(genome.at(x + kmersize - i - 1));
			} else { //!detx
				if (genome.at(x+i) == revcomp(genome.at(x+kmersize -i - 1))) {
					xchar = genome.at(x+i);
				} else {
					detx = true;
					flipx = (genome.at(x+i) > revcomp(genome.at(x+kmersize -i - 1)));
					if (flipx) {
						xchar = revcomp(genome.at(x+kmersize -i - 1));
					} else {
						xchar = genome.at(x + i);
					}
				}
			}
			if (dety && !flipy) {
				ychar = genome.at(y+i);
			} else if (dety && flipy) {
				ychar = revcomp(genome.at(y + kmersize - i - 1));
			} else { //!dety
				if (genome.at(y+i) == revcomp(genome.at(y+kmersize -i - 1))) {
					ychar = genome.at(y+i);
				} else {
					dety = true;
					flipy = (genome.at(y+i) > revcomp(genome.at(y+kmersize -i - 1)));
					if (flipy) {
						ychar = revcomp(genome.at(y+kmersize -i - 1));
					} else {
						ychar = genome.at(y + i);
					}
				}
			}

			if (xchar != ychar) return xchar < ychar;
		}
		return false;
	}

	string genome;
	long genSize;
	int kmersize;
};

void usage(int argc, char * argv[]) {
	exit(1);
}

int main(int argc, char * argv[]) {

	char ch;
	string genomeBase;
	string genome;
	int kmersize = -1;
	while ((ch = getopt(argc, argv, "k:g:")) != -1) {
		switch (ch) {
			case 'g':
				genomeBase = optarg;
				genome = read_genome(genomeBase + ".fa");
				break;
			case 'k':
				kmersize = atoi(optarg);
				break;
		}
	}
	//string task = argv[optind];
	if (genome == "" || kmersize == -1) { 
		cerr << "Invalid Params.\n";
		usage(argc, argv);
	}
	suffixArray sa(genome);
	sa.load(genomeBase + ".sa");


	//deque<string> badKmers; //kmers that don't hit ecoli
	deque<long> badKmers; //kmers that don't hit ecoli
	deque<string> goodKmers; //kmers that hit ecoli
	deque<string> genomicKmers; //kmers in the genome
	vector<bool> hits(genome.size(), false);

	string readsPseudoGenome;
	string seq;
	cout << "Mapping kmers...\n";
	while (getline(cin, seq)) {
		if (seq.length() < kmersize) continue;
		int hit = max(sa.find(seq), sa.find(revcomp(seq)));
		if (hit != -1) {
			for (int i = 0; i < seq.length() - kmersize + 1; i++) {
				hits.at((hit + i) % genome.size()) = true;
			}
		} else {
			for (int i = 0; i < seq.length() - kmersize + 1; i++) {
				string kmer = seq.substr(i,kmersize);
				int hit = max(sa.find(kmer), sa.find(revcomp(kmer)));
				if (hit != -1) {
					hits.at(hit) = true;
				} else {
					//badKmers.push_back(rcnorm(kmer));
					//if (seq.find_first_of("N", i, kmersize) == string::npos) {
					badKmers.push_back(readsPseudoGenome.size() + i);
					//}
				}
			}
			readsPseudoGenome += seq + '$';
		}
	}
	for (int i = 0; i < genome.size(); i++) {
		string kmer = genome.substr(i,kmersize);
		genomicKmers.push_back(rcnorm(kmer));
		if (hits.at(i)){
			goodKmers.push_back(rcnorm(kmer));
		}
	}
	cout << "Sorting goodKmers (" << goodKmers.size() << " elements)...\n";
	sort(goodKmers.begin(), goodKmers.end());
	goodKmers.resize(unique(goodKmers.begin(), goodKmers.end()) - goodKmers.begin());
	cout << "Number of distinct genomic kmers in the data: " << goodKmers.size() << endl;
	cout << "Sorting badKmers (" << badKmers.size() << " elements)...\n";
	/*ofstream dump;
	open_file(dump, "badKmers.txt");
	for (long i = 0; i < badKmers.size(); i++) {
		dump << badKmers.at(i) << endl;
	}
	dump.close();
	*/
	compltNum2NumFunctor complt;
	complt.genome = readsPseudoGenome;
	complt.kmersize = kmersize;
	complt.genSize = readsPseudoGenome.size();
	sort(badKmers.begin(), badKmers.end(), complt);
	compeqNum2NumFunctor compeq;
	compeq.genome = readsPseudoGenome;
	compeq.kmersize = kmersize;
	badKmers.resize(unique(badKmers.begin(), badKmers.end(), compeq) - badKmers.begin());

	/*sort(badKmers.begin(), badKmers.end());
	  badKmers.resize(unique(badKmers.begin(), badKmers.end()) - badKmers.begin());
	 */
	cout << "Number of distinct kmers in the data: " << goodKmers.size() + badKmers.size() << endl;
	double spec = double(goodKmers.size()) / double(goodKmers.size() + badKmers.size());
	cout << "Specificity: " << spec << endl;
	cout << "Sorting genomicKmers (" << genomicKmers.size() << " elements)...\n";
	sort(genomicKmers.begin(), genomicKmers.end());
	genomicKmers.resize(unique(genomicKmers.begin(), genomicKmers.end()) - genomicKmers.begin());
	cout << "Number of distinct kmers in the genome: " << genomicKmers.size() << endl;
	double sens = double(goodKmers.size()) / double(genomicKmers.size());
	cout << "Sensitivity: " << sens << endl;
	cout << "RES\t" << goodKmers.size() + badKmers.size() << "\t" << goodKmers.size() << "\t" << spec << "\t" << sens << endl;


}

