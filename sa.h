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
#include "defs.h"

struct sah__compltNum2NumFunctor {
	bool operator()(int x, int y) const {
		if (x==y) return false;
		for (int i = 0; i < genSize - max(x,y) + 1; i++) {
			if (genome[x+i] != genome[y+i]) return genome[x+i] < genome[y+i];
		}
		return x > y; //no, this is not a mistake
	}
	string genome;
	int genSize;
};

struct sah__compltNum2StrFunctor {
	bool operator()(int x, string y) const {
		return genome.substr(x, y.length()) < y;
	}
	string genome;
	int genSize;
};




class suffixArray {
	public:
		int genSize;
		vector<int> sa;

		suffixArray(string g) {
			genSize = g.size();
			genome = g + g;
		}
		void build() {
			cerr << "Building suffix array...\n";
			for (int i = 0; i < genSize; i++) {
				sa.push_back(i);
			}

			sah__compltNum2NumFunctor complt; 
			complt.genome = genome;
			complt.genSize = genSize;
			sort(sa.begin(), sa.end(), complt);
		}
		/*void build(int kmersize, char skip) {
			cerr << "Building suffix array...\n";
			for (int i = 0; i < genSize; i++) {
				if (genome.find_first_of('$', i, kmersize) != string::npos) continue;
				sa.push_back(i);
			}

			sah__compltNum2NumFunctor complt; 
			complt.genome = genome;
			complt.genSize = genSize;
			sort(sa.begin(), sa.end(), complt);
		}
		*/

		void save(string filename, int k = 0) {
			ofstream out;
			open_file(out, filename);
			if (k==0) {
				for (int i = 0; i < sa.size(); i++) out << sa[i] << endl;
			} else {
				for (int i = 0; i < sa.size(); i++) out << sa[i] << "\t" << genome.substr(sa[i],k) <<  endl;
			}
			out.close();
		}
		void load(string filename) {
			sa.clear();
			ifstream in;
			open_file(in, filename);
			string sbuf;
			while (getline(in, sbuf)) {
				istringstream line(sbuf);
				int val;
				line >> val;
				sa.push_back(val);
			}
			in.close();
		}
		int find(string kmer) {
			sah__compltNum2StrFunctor complt; 
			complt.genome = genome;
			complt.genSize = genSize;
			vector<int>::iterator it = lower_bound(sa.begin(), sa.end(), kmer, complt);
			if (it == sa.end() || (genome.substr(*it, kmer.length()) != kmer) ) {
				return -1;
			}
			return *it;
		}
		void find(string kmer, vector<int> &hits) {
			hits.clear();
			sah__compltNum2StrFunctor complt; 
			complt.genome = genome;
			complt.genSize = genSize;
			vector<int>::iterator it = lower_bound(sa.begin(), sa.end(), kmer, complt);
			while (it != sa.end() || (genome.substr(*it, kmer.length()) == kmer) ) {
				hits.push_back(*it);
				it++;
			}
			return;
		}

private:
	string genome;

};

