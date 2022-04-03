#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

using namespace std;
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream> 
#include <time.h>
#include <sys/time.h>

const double EPSILON = 0.00001; // small constant
const int M = 1000000;			// big constant

struct Instance
{
	int m; 						// number of bar chart types
	int c; 						// capacity of the bins
	int UB;						// quick upper bound
	vector<vector<int> > BCs;	// bar charts
	void print();
};

struct Solution
{
	int opt, LB, UB, Nvar, Nconstr, Ncoeff, Ncuts;
	double timeP, timeT;
    vector<vector<int>> assignedBin;
};

double getCPUTime();
bool sortBC(const vector<int>& v1, const vector<int>& v2);
Instance readInstance(string filename);
void printInfo(const string& pathAndFileout, const Solution& sol, const string& filein);
void printSInfo(const Solution& sol, const Instance& inst);

#endif 