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

class Line {
public:
	string idx;
	string h1;
	string h2;
};

int types[6] = {0,0,0,0,0,0};
int bigblock = 0;


void process_block(deque<Line> & block) {
	int type;
	int cats[4] = {0,0,0,0}; 
	// 0 = -1 -1, 
	// 1 = -1 maps, 
	// 2 = maps -1, 
	// 3 = maps maps
	for (int i = 0; i < block.size(); i++) {
		if (block[i].h1 == "-1" && block[i].h2 == "-1") {
			cats[0]++;
		} else if (block[i].h1 == "-1" && block[i].h2 == "maps") {
			cats[1]++;
		} else if (block[i].h1 == "maps" && block[i].h2 == "-1") {
			cats[2]++;
		} else if (block[i].h1 == "maps" && block[i].h2 == "maps") {
			cats[3]++;
		}
	}
	if (cats[1] + cats[3] > 0) {
		if (cats[3] > 1) {
			type = 0;
		} else if (cats[3] == 1) {
			type = 1;
		} else if (cats[3] == 0) {
			type = 2;
		} 
	} else {
		if (cats[2] > 0) {
			type = 3;
		} else {
			type = 4;
		}
	}


	//for (int i = 0; i < block.size(); i++) cout << "\t" << block[i].idx << "\t" << block[i].h1 << "\t" << block[i].h2 << endl; cout << "CAT\t" <<type << endl <<  endl;

	types[type]++;
}


int main(int argc, char * argv[]) {
	deque<Line> block;
	Line line, lastline;
	cin >> line.idx >> line.h1 >> line.h2;
	lastline = line; 
	block.push_back(line);
	while (cin >> line.idx >> line.h1 >> line.h2) {
		if (lastline.idx == line.idx) {
			block.push_back(line);
		} else {
			process_block(block);
			block.clear();
			block.push_back(line); 
		}
		lastline = line; 
	}
	process_block(block); 
	cout << "0\t" << types[0] << endl;
	cout << "1\t" << types[1] << endl;
	cout << "2\t" << types[2] << endl;
	cout << "3\t" << types[3] << endl;
	cout << "4\t" << types[4] << endl;
	cout << "unknown\t" << types[5] << endl;
	cout << "bigblock\t" << bigblock << endl;

	return 0;
}

/*
   int main(int argc, char * argv[]) {
   deque<Line> block;
   Line line, lastline;
   cin >> line.count >> line.idx >> line.h1 >> line.h2;
   lastline = line; 
   block.push_back(line);
   while (cin >> line.count >> line.idx >> line.h1 >> line.h2) {
   if (lastline.idx == line.idx) {
   block.push_back(line);
   } else {
   process_block(block);
   block.clear();
   block.push_back(line); 
   }
   lastline = line; 
   }
   process_block(block); 
   cout << "0\t" << types[0] << endl;
   cout << "1\t" << types[1] << endl;
   cout << "2\t" << types[2] << endl;
   cout << "3\t" << types[3] << endl;
   cout << "4\t" << types[4] << endl;
   cout << "bigblcok\t" << bigblock << endl;

   return 0;
   }

   void process_block(deque<Line> & block) {
   if (block.size() > 2) {
   bigblock++;
   return;
   }


   int type;
   if (block[0].h2 == "-1" ) { //type3 or 4
   if (block.size() == 1) {
   type = 4;
   } else {
   type = 3;
   }
   } else if (block.size() == 1) {
   if (block[0].h1 == "-1") {
   type = 2;
   } else {
   type = 0; 
   }
   } else if (block[1].count == 1) {
   type = 1; 
   } else {
   type = 0;
   }
//for (int i = 0; i < block.size(); i++) cout << block[i].count << "\t" << block[i].idx << "\t" << block[i].h1 << "\t" << block[i].h2 << endl; cout << "CAT\t" <<type << endl <<  endl;

types[type]++;
}



*/
