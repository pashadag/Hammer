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


int main(int argc, char * argv[]) {

	string header, seq, junk, qv;
	while (getline(cin, header)) {
		getline(cin, seq);
		getline(cin, junk);
		getline(cin, qv);
		if (seq.find('N') == string::npos)  cout << header << endl << seq << endl << junk << endl << qv << endl;
	}
}


