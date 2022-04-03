#include "helper_functions.h"

void Instance::print(){
	cout << "m = " << m << " " << "c = " << c << " " << "UB = " << UB << endl;
	for(int i=0; i<m;i++)
		cout << BCs[i][0] << " " << BCs[i][1] << " " << BCs[i][2]  << endl;
}

double getCPUTime(){
	return (double)clock() / CLOCKS_PER_SEC;
}

bool sortBC(const vector<int>& v1, const vector<int>& v2) {
    return v1[0] +  v1[1] >  v2[0] + v2[1];
}

Instance readInstance(string filename) {
    // define variables
	Instance inst;

    // open the file
    ifstream file(filename); 

    // read the file
    if (file.is_open()) { //if the file is open
        string line;

        // first line contains number of BCs
        getline(file, line, '\n');
        inst.m = stoi(line);

        // reshape the array that will hold the bar heights
        inst.BCs.resize(inst.m, vector<int>(2, 0));

        // second line contains total cell height
        getline(file, line, '\n');
        inst.c = stoi(line);

        // the remaining lines contain two bar heights each
        for (int i = 0; i < inst.m; i++) {
            getline(file, line, ',');
            inst.BCs[i][0] = stoi(line);
            getline(file, line, ',');
            inst.BCs[i][1] = stoi(line);
            getline(file, line, '\n');
            inst.BCs[i][2] = stoi(line);
        }

        // close the file
        file.close(); 
    }
    else {
        // if the file cannot be opened: print error and return default Inst
        cout << "Unable to open file"; 
    }

	sort(inst.BCs.begin(), inst.BCs.end(), sortBC);
	
	// First Fit decreasing
	vector<int> filling (M,0);
	for(int i = 0; i <inst.m;i++){
		int idx = 0;
		for(int j = 0; j <inst.BCs[i][2];j++){
			while (filling[idx] + inst.BCs[i][0] > inst.c || filling[idx +1] + inst.BCs[i][1] > inst.c) idx++;
			filling[idx] += inst.BCs[i][0]; filling[idx+1] += inst.BCs[i][1];
		}
	}
	inst.UB = 0; while(filling[inst.UB] > 0) inst.UB++;
	
    return inst;
}

void printInfo(const string& pathAndFileout, const Solution& sol, const string& filein){
	string nameFile = pathAndFileout;
	std::ofstream file(nameFile.c_str(), std::ios::out | std::ios::app);
	file << filein << "\t"  << sol.opt << "\t" << sol.timeT << "\t" << sol.timeP << "\t" << sol.LB << "\t" << sol.UB  << "\t" << sol.Nvar << "\t" << sol.Nconstr << "\t" << sol.Ncoeff  << "\t" << sol.Ncuts << endl;
	file.close();
}

void printSInfo(const Solution& sol, const Instance& inst){
	vector<vector<int> > bins1(sol.UB);
	vector<vector<int> > bins2(sol.UB);
	vector<int>	load(sol.UB,0);
	for (int i = 0; i < sol.assignedBin.size();i++){
		for (int j = 0; j < sol.assignedBin[i].size();j++){
			load[sol.assignedBin[i][j]] += inst.BCs[i][0];
			load[sol.assignedBin[i][j]+1] += inst.BCs[i][1];
			bins1[sol.assignedBin[i][j]].push_back(i);
			bins2[sol.assignedBin[i][j]+1].push_back(i);
		}
	}
	for (int i = 0; i < sol.UB;i++){
		cout << "Bin " << i << "\t" << load[i] << "\t [";
		for (int j = 0; j <  bins1[i].size(); j++) cout << bins1[i][j] << " ";
		cout << "] [";
		for (int j = 0; j <  bins2[i].size(); j++) cout << bins2[i][j] << " ";
		cout << "]" << endl;
	}
}