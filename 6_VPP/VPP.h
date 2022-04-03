#ifndef NAIVE_H
#define NAIVE_H

using namespace std;
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <math.h> 
#include "gurobi_c++.h"
#include "helper_functions.h"

struct Graph
{
    vector<vector<int>> A;
    int Narcs, Nvert;
	vector<vector<int>> ITN;
	vector<vector<int>> NTI;
};

Solution VPP(const Instance& inst);
Graph makeGraphVP(const Instance& inst);

#endif 