#include "GALO.h"

Solution GALO(const Instance& inst) {

	// local variable
	double start = getCPUTime(); 			// starting time
    Solution sol; 							// solution to return
	sol.assignedBin.resize(inst.m);
	
	// GALO
	vector<int> filling (inst.UB*2,0);
	vector<int> remaining;
	int nbP = 0; int mini = 0;
	for(int i = 0; i <inst.m;i++) {
		nbP += inst.BCs[i][2];
		remaining.push_back(inst.BCs[i][2]);
	}
	
	while(nbP > 0){
		bool hbf = false;
		for(int j = 0; j<inst.m;j++){
			if(remaining[j] == 0 || filling[mini] + inst.BCs[j][0] > inst.c || filling[mini +1] + inst.BCs[j][1] > inst.c) 
				continue;
			filling[mini] += inst.BCs[j][0]; filling[mini+1] += inst.BCs[j][1]; remaining[j]--; nbP--; sol.assignedBin[j].push_back(mini);
			hbf = true;
			break;
		}
		if (!hbf) mini++;
	}
	sol.UB = 0; while(filling[sol.UB] > 0) sol.UB++;
	sol.timeT = getCPUTime() - start;
	 
	sol.timeP = -1;
    sol.Nvar = -1;
    sol.Nconstr = -1;
    sol.Ncoeff = -1;
	sol.opt = -1;
	sol.LB = -1;
	
    return sol;
}

