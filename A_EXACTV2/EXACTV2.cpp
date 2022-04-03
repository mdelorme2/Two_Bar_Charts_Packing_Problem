#include "EXACTV2.h"

Graph makeGraph2(const Instance& inst) {
	
	// Local variables
	Graph G;
    
	// Initialisation  
	G.NTI.resize(inst.c + 1, vector<int>(inst.c + 1, -1));
 
    int newC1, newC2, newC3, newC4;
	 
	// First node (0,0) is active
    G.NTI[0][0] = 0;               		
	G.ITN.push_back({0,0}); 
    G.Nvert = 1;                   		

	// Loop over all BCs
    for (int k = 0; k < inst.m; k++) {  	
 		int nbNodes = G.Nvert;              			// fix the number of nodes for that iteration
		vector<bool> alreadyVisited (nbNodes,false);
		for (int i = 0; i < nbNodes; i++) {     		// loop over all the existing vertices
			for (int j = 0; j < inst.BCs[k][2]; j++) { 	// for each unit of demand of BC type k
				newC1 = G.ITN[i][0] + j * inst.BCs[k][0];	
				newC2 = G.ITN[i][1] + j * inst.BCs[k][1];
				newC3 = G.ITN[i][0] + (j+1) * inst.BCs[k][0];	
				newC4 = G.ITN[i][1] + (j+1) * inst.BCs[k][1];
				int idxN = G.NTI[newC1][newC2];
				if(idxN >= 0 && idxN < nbNodes){
					if (alreadyVisited[idxN] == true) break;
					else alreadyVisited[idxN] = true;
				}
				if(newC3 <= inst.c && newC4 <= inst.c){ // if one unit of the BC type can be included
					if (G.NTI[newC3][newC4] == -1) {	// if the new head does not exist yet, create it
						G.NTI[newC3][newC4] = G.Nvert;	
						G.ITN.push_back({newC3,newC4});
						G.Nvert++; 			
					}
					G.A.push_back({G.NTI[newC1][newC2],G.NTI[newC3][newC4],k});   
				}
				else break;
			}
		}
    }
	
	// Add transition arcs
    for (int i = 0; i < G.Nvert; i++) G.A.push_back({i, 0, -1}); 
	
	// Erase duplicates in arcs
	set<vector<int>> s(G.A.begin(),G.A.end());
	G.A.assign(s.begin(),s.end());
	G.Narcs = G.A.size();
	cout << "---------------------------------------" << endl;
	
    // return the Graph object
    return G;
}

mycallbackT2::mycallbackT2(const Graph& xG, const vector<GRBVar>& xfa, const vector<GRBVar>& xla, const vector<vector<int>>& xlossArcs, const int& xUB) {
    // initialize the Callback object
    this->Ncuts = 0;  // the number of cuts, initially 0
    this->G = xG;     // the arcflow graph
    this->fa = xfa;   // the flow
    this->la = xla;   // the loss arcs
	this->lossArcs = xlossArcs;	  
	this->UB = xUB;
}

int mycallbackT2::getNcuts() {
    return this->Ncuts; // return the number of cuts
}

void mycallbackT2::callback() {
    if (where == GRB_CB_MIPSOL) { // if a new incumbent is found
	
        // find the incumbent solution
        vector<int> f(G.Narcs); vector<int> l(lossArcs.size());                 			
        for (int i = 0; i < G.Narcs; i++) f[i] = ceil(getSolution(fa[i]) - EPSILON); 
		for (int i = 0; i < lossArcs.size(); i++) l[i] = ceil(getSolution(la[i]) - EPSILON); 
		
        // find all cycles
		int c = lossArcs.back()[1];
		vector<GRBLinExpr> tOut (c+1, 0);   								// the amount of flow leaving each transition vertex
		vector<vector<GRBLinExpr>> tArcs (c+1,vector<GRBLinExpr>(c+1,0));   // the arcs associated with each metaArc
		vector<vector<vector<int>>> tUsed (c + 1);
		vector<vector<int>> arcsUsed (G.Nvert);
		int nbRem = 0;
		for (int i = 0; i < G.Narcs; i++) {
			if (G.A[i][2] == -1){
				tOut[c - G.ITN[G.A[i][0]][0]] += fa[i];
				tArcs[c - G.ITN[G.A[i][0]][0]][G.ITN[G.A[i][0]][1]] += fa[i];
			}
			for (int k = 0; k < f[i]; k++){  
				if(G.A[i][2] >= 0) arcsUsed[G.A[i][0]].push_back(i);
				else{
					arcsUsed[G.A[i][0]].push_back(i);
					tUsed[c - G.ITN[G.A[i][0]][0]].push_back({c - G.ITN[G.A[i][0]][0],G.ITN[G.A[i][0]][1],0});
					nbRem++;
				}
			}
		}
		for (int i = 0; i < lossArcs.size(); i++) { 
			tOut[lossArcs[i][0]] += la[i];
			tArcs[lossArcs[i][0]][lossArcs[i][1]] += la[i];
			for (int k = 0; k < l[i]; k++){  
				tUsed[lossArcs[i][0]].push_back({lossArcs[i][0],lossArcs[i][1],-1});
				nbRem++;
			}
		}
		
		// get the metaCycles
		vector<vector<vector<int>>> metaCycles;
		while(nbRem > 0){
			vector<vector<int>> metaCycle; int tail = 0;
			while(tUsed[tail].size() == 0) tail++;
			// Create cycle
			while(tUsed[tail].size() > 0){
				metaCycle.push_back(tUsed[tail].back()); 
				tUsed[tail].pop_back(); 
				tail = metaCycle.back()[1];
				nbRem--;
			}
			// Extend cycle
			for(int i = 0; i < metaCycle.size();i++){
				if(tUsed[metaCycle[i][0]].size() > 0){
					tail = metaCycle[i][0];
					int tempI = i;
					while(tUsed[tail].size() > 0){
						metaCycle.insert(metaCycle.begin() + tempI,tUsed[tail].back()); tempI++;
						tUsed[tail].pop_back(); 
						tail = metaCycle[tempI-1][1];
						nbRem--;
					}
				}
			}
			metaCycles.push_back(metaCycle);
		}

		for(int i = 0; i < metaCycles.size();i++){
			cout << "MetaCycle " << i << ": "; 
			for(int j = 0; j < metaCycles[i].size();j++){
				cout << "( ";
				for(int k = 0; k < metaCycles[i][j].size();k++){
					cout << metaCycles[i][j][k] << " ";
				}
				cout << ") ";
			}
			cout << endl;
		}
		
        // for each supplementary metaCycle, add a no-good cut
		for(int i = 1; i < metaCycles.size();i++){	
			int count = 0;
			GRBLinExpr LHS = 0;
			vector<bool> isNodeInTheCut (c+1,false);	
			for(int j = 0; j < metaCycles[i].size();j++) isNodeInTheCut[metaCycles[i][j][0]] = true; 
			GRBLinExpr RHS = 0;
			for (int k = 0; k < c+1; k++) { 
				for (int l = 0; l < c+1; l++) { 
					if(isNodeInTheCut[k] && !isNodeInTheCut[l]){
						RHS += tArcs[k][l];
					}
					if(isNodeInTheCut[k] && isNodeInTheCut[l]){
						LHS += tArcs[k][l]; 
						count++;
					}
				}
			}
			// add the cut to the set of variables
			addLazy(LHS <= RHS * UB * count);  			
			Ncuts++;     
		}
		cout << "Added " << Ncuts << " cuts " << endl;
		cout << "-------------------------------------------------------" << endl;
    }
}

Solution EXACTV2(const Instance& inst) {

	// local variable
	double start = getCPUTime(); 				// starting time
	Solution sol; sol.Ncuts = 0;				// solution to return
	Graph G = makeGraph2(inst);   
	vector<bool> activeNodes (inst.c+1,false);	// for loss arcs
	
    // create a model
    GRBEnv env = GRBEnv();              	// create an environment
    GRBModel model = GRBModel(env);         // create a new model

    // declaration of the variables for the model
    vector<GRBVar> fa(G.Narcs);
	
    // initizalization of the variables for the model
    for (int i = 0; i < G.Narcs; i++) { 
		activeNodes[inst.c - G.ITN[G.A[i][0]][0]] = true; activeNodes[G.ITN[G.A[i][0]][1]] = true;
		//cout << G.ITN[G.A[i][0]][0] << "," << G.ITN[G.A[i][0]][1] << " " << G.ITN[G.A[i][1]][0] << "," << G.ITN[G.A[i][1]][1] << " " << G.A[i][2] << endl;
		if(G.A[i][2] != -1) fa[i] = model.addVar(0, inst.BCs[G.A[i][2]][2], 0, GRB_INTEGER);
		else fa[i] = model.addVar(0, inst.UB, 0, GRB_INTEGER);
    }
	
	// create loss arcs structure
	vector<vector<int>> lossArcs;
	activeNodes[inst.c] = true;
	vector<int> isActive; for(int i = 0; i < inst.c+1;i++) if (activeNodes[i]) isActive.push_back(i);
	for (int i = 0; i < isActive.size()-1;i++) lossArcs.push_back({isActive[i],isActive[i+1]});
	vector<GRBVar> la(lossArcs.size());
	for (int i = 0; i < lossArcs.size(); i++){
		// cout << lossArcs[i][0] << " " << lossArcs[i][1] << endl;
		la[i] = model.addVar(0, inst.UB, 0, GRB_INTEGER);
	}
    model.update();	

	// declare linear expressions
    vector<GRBLinExpr> fIn (G.Nvert, 0);    		// the amount of flow entering each vertex
    vector<GRBLinExpr> fOut (G.Nvert, 0);   		// the amount of flow leaving each vertex
    vector<GRBLinExpr> tIn (inst.c+1, 0);    		// the amount of flow entering each transition vertex
    vector<GRBLinExpr> tOut (inst.c+1, 0);   		// the amount of flow leaving each transition vertex
    vector<GRBLinExpr> typeUsed(inst.m + 1, 0); 	// the amount of arcs used of each BC type 
    
	// calculate the linear expressions	
    for (int i = 0; i < G.Narcs; i++) {   		// loop over all arcs
        fIn[G.A[i][1]] += fa[i];        		// inflow
        fOut[G.A[i][0]] += fa[i];       		// outflow
		if (G.A[i][2] == -1){
			tOut[inst.c - G.ITN[G.A[i][0]][0]] += fa[i];
			tIn[G.ITN[G.A[i][0]][1]] += fa[i];
		}
        typeUsed[G.A[i][2]+1] += fa[i];    		// number of BCs used of certain type
    }
	for (int i = 0; i < lossArcs.size(); i++){
		tOut[lossArcs[i][0]] += la[i];
		tIn[lossArcs[i][1]] += la[i];
	}
	
	model.update();

    // create flow conservation constraints 
    for (int v = 0; v < G.Nvert; v++) {       	// loop over all vertices
        model.addConstr(fIn[v] == fOut[v]); 	// inflow = outflow
    }
    for (int v = 0; v < isActive.size(); v++) { // loop over all transition vertices	
		model.addConstr(tIn[isActive[v]] == tOut[isActive[v]]);
	}

    // create constraint ensuring at least one cycle containing (0, 0)
    model.addConstr(fOut[0] >= 1); model.addConstr(tOut[0] >= 1); 

    // create BC type constraints
    for (int i = 0; i < inst.m; i++) {           // loop over all BC types
		if (inst.BCs[i][2] == 1) model.addConstr(typeUsed[i + 1] == inst.BCs[i][2]); 	// demand met exactly if 1
		else model.addConstr(typeUsed[i + 1] >= inst.BCs[i][2]); 						// demand met otherwise
    }

    // set the objective: minimize the number of loss arcs used
    model.setObjective(typeUsed[0], GRB_MINIMIZE);

    // change some settings
    model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
    model.getEnv().set(GRB_IntParam_Threads, 1);
    model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);

    // find the optimal solution
	sol.timeP = getCPUTime() - start;			// preprocessing time

    // create a callback
    mycallbackT2 cb = mycallbackT2(G, fa, la, lossArcs,inst.UB);
    model.set(GRB_IntParam_LazyConstraints, 1); 		// indicate that we want to add Lazy Constraints
    model.setCallback(&cb);                     		// link the callback to the model
	
    model.optimize();

    // store the results in a Solution object
    sol.Nvar = model.get(GRB_IntAttr_NumVars);       // number of variables
    sol.Nconstr = model.get(GRB_IntAttr_NumConstrs); // number of constraints
    sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);      // number of non-zero coefficients		
	sol.opt = 0;
	sol.LB = ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
	sol.UB = -1;
	sol.Ncuts = cb.getNcuts();
	
    if (model.get(GRB_IntAttr_SolCount) >= 1) { 		// if a solution has been found
         sol.UB = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		
		// get bin for each bar chart
		vector<vector<vector<int>>> tUsed (inst.c + 1);
		vector<vector<int>> arcsUsed (G.Nvert);
		int nbRem = 0;
		for (int i = 0; i < G.Narcs; i++) {
			for (int k = 0; k < ceil(fa[i].get(GRB_DoubleAttr_X) - EPSILON); k++){  
				if(G.A[i][2] >= 0) arcsUsed[G.A[i][0]].push_back(i);
				else{
					arcsUsed[G.A[i][0]].push_back(i);
					tUsed[inst.c - G.ITN[G.A[i][0]][0]].push_back({inst.c - G.ITN[G.A[i][0]][0],G.ITN[G.A[i][0]][1],0});
					nbRem++;
				}
			}
		}
		for (int i = 0; i < lossArcs.size(); i++) { 
			for (int k = 0; k < ceil(la[i].get(GRB_DoubleAttr_X) - EPSILON); k++){  
				tUsed[lossArcs[i][0]].push_back({lossArcs[i][0],lossArcs[i][1],-1});
				nbRem++;
			}
		}
		
		// get the metaCycles
		vector<vector<vector<int>>> metaCycles;
		while(nbRem > 0){
			vector<vector<int>> metaCycle; int tail = 0;
			while(tUsed[tail].size() == 0) tail++;
			// Create cycle
			while(tUsed[tail].size() > 0){
				metaCycle.push_back(tUsed[tail].back()); 
				tUsed[tail].pop_back(); 
				tail = metaCycle.back()[1];
				nbRem--;
			}
			// Extend cycle
			for(int i = 0; i < metaCycle.size();i++){
				if(tUsed[metaCycle[i][0]].size() > 0){
					tail = metaCycle[i][0];
					int tempI = i;
					while(tUsed[tail].size() > 0){
						metaCycle.insert(metaCycle.begin() + tempI,tUsed[tail].back()); tempI++;
						tUsed[tail].pop_back(); 
						tail = metaCycle[tempI-1][1];
						nbRem--;
					}
				}
			}
			metaCycles.push_back(metaCycle);
		}
		
		for(int i = 0; i < metaCycles.size();i++){
			cout << "MetaCycle " << i << ": "; 
			for(int j = 0; j < metaCycles[i].size();j++){
				cout << "( ";
				for(int k = 0; k < metaCycles[i][j].size();k++){
					cout << metaCycles[i][j][k] << " ";
				}
				cout << ") ";
			}
			cout << endl;
		}
		
		// find paths
		vector<vector<vector<vector<int>>>> paths(inst.c+1); 
		for(int i =0; i<inst.c+1;i++) paths[i].resize(inst.c+1);
		while(arcsUsed[0].size() > 0){
			vector<int> path; int tail = 0;
			while(tail != 0 || path.size() == 0){
				path.push_back(arcsUsed[tail].back()); 
				arcsUsed[tail].pop_back(); 
				tail = G.A[path.back()][1];
			}
			paths[inst.c - G.ITN[G.A[path.back()][0]][0]][G.ITN[G.A[path.back()][0]][1]].push_back(path);
		}

		if(metaCycles.size() == 1){
			cout << "No subtour, solution is optimal" << endl;
			// get bin for each BC type
			sol.assignedBin.resize(inst.m);
			int binIdx = 0;
			for(int i = 0; i < metaCycles[0].size();i++){
				if(metaCycles[0][i][2] == -1) continue;
				vector<int> path = paths[metaCycles[0][i][0]][metaCycles[0][i][1]].back();  paths[metaCycles[0][i][0]][metaCycles[0][i][1]].pop_back();
				for (int j = 0; j< path.size();j++) if(G.A[path[j]][2] >= 0) sol.assignedBin[G.A[path[j]][2]].push_back(binIdx);
				binIdx++;
			}
			if(sol.LB == sol.UB) sol.opt = 1;
		}
    }

    sol.timeT = getCPUTime() - start;	
    return sol;
}
