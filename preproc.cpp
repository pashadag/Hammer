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
#include<cassert>
#include "defs.h"
using namespace std;

int qvoffset;
string mode;

double oct2phred(string qoct)  {
	float freq = 1;	
	for (int i = 0; i < qoct.length(); i++) {
		freq *= 1 - pow(10, -float(qoct[i] - qvoffset)/10.0);
	}
	return freq;
}               

void trimbs(string & seq, string & qv, int starttrim, int endtrim) {
	assert(seq.length() == qv.length());
	char lowqual = qvoffset + 2; // B or #
	if (endtrim == -1) endtrim = qv.length();
	endtrim = min(endtrim, (int) qv.length());
	for (int i = starttrim; i < endtrim; i++) {
		if (qv[i] != lowqual) break; else starttrim++; 
	}
	for (int i = endtrim - 1; i >= 0;  i--) {
		if (qv[i] != lowqual) break; else endtrim--;
	}
	if (endtrim > starttrim) {
		seq = seq.substr(starttrim, endtrim - starttrim  );
		qv = qv.substr(starttrim, endtrim - starttrim  );
	} else {
		seq = ""; 
		qv  = "";
	}
	return;
}

int main(int argc, char * argv[]) {
	assert(argc == 5);
	int starttrim = atoi(argv[1]); 
	int endtrim = atoi(argv[2]);  //number of characters to keep
	qvoffset = atoi(argv[3]);
	mode = argv[4];
	if (mode != "raw" && mode != "fastq" ) {
		cerr << "Invalid mode : " << mode << endl;
		exit(1);
	}

	string header, seq, junk, qv;
	while (getline(cin, header)) {
		getline(cin, seq);
		getline(cin, junk);
		getline(cin, qv);
		trimbs(seq, qv, starttrim, endtrim);
		if (seq != "") {
			if (mode == "raw" ) {
				cout << seq << "\t" << oct2phred(qv)<< endl;
			} else if (mode == "fastq") {
				cout << header << endl << seq << endl << junk << endl << qv << endl;
			}
		}
	}
}


