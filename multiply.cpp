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
ostringstream oMsg;



int main(int argc, char * argv[]) {
	int maxBound = 1000000;
	if (argc > 1) maxBound = atoi(argv[1]);
	string kmer;
	int count;
	while (cin >> kmer) {
		cin >> count;
		if (count > maxBound) count = maxBound;
		for (int i = 0; i < count; i++) {
			cout << ">read\n" << kmer << endl;
		}
	}
}

