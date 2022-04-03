#include "EXACTV1.h"

Graph makeGraph1(const Instance& inst) {

	// Local variables
	Graph G;
	vector<bool> Sbool(inst.c + 1, false); // set of reachable start vertices S
    
	// Initialisation  
	G.NTI.resize(inst.c + 1, vector<int>(inst.c + 1, -1));
    Sbool[0] = true;
 
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
						Sbool[newC4] = true;
					}
					G.A.push_back({G.NTI[newC1][newC2],G.NTI[newC3][newC4],k});   
				}
				else break;
			}
		}
    }

    // Find all reachable vertices by 'shifting' the base layer
	int nbArcs = G.A.size();							// fix the number of arcs for all iteration
	set<vector<int>> se(G.A.begin(),G.A.end());
    for (int s = 1; s < inst.c+1; s++){
		if (Sbool[s]){
			for(int a =0; a<nbArcs;a++){
				if (G.ITN[G.A[a][1]][0] + s <= inst.c) {
					newC1 = G.ITN[G.A[a][0]][0] + s;	
					newC2 = G.ITN[G.A[a][0]][1];
					if (G.NTI[newC1][newC2] == -1) {	// if the new tail does not exist yet, create it
						G.NTI[newC1][newC2] = G.Nvert;	
						G.ITN.push_back({newC1,newC2});
						G.Nvert++; 		
					}
					newC3 = G.ITN[G.A[a][1]][0] + s;	
					newC4 = G.ITN[G.A[a][1]][1];
					if (G.NTI[newC3][newC4] == -1) {	// if the new head does not exist yet, create it
						G.NTI[newC3][newC4] = G.Nvert;	
						G.ITN.push_back({newC3,newC4});
						G.Nvert++; 		
					}
					se.insert({G.NTI[newC1][newC2],G.NTI[newC3][newC4],G.A[a][2]});   
				}
			}
		}
    }

    // Add the loss arcs 
    for (int i = 1; i < G.Nvert; i++) {     	
		newC1 = G.ITN[i][1];	
		newC2 = 0;
		if (G.NTI[newC1][newC2] == -1) {	// if the new head does not exist yet, create it
			G.NTI[newC1][newC2] = G.Nvert;	
			G.ITN.push_back({newC1,newC2});
			G.Nvert++; 		
		}
		se.insert({i, G.NTI[newC1][newC2], -1}); 
    }
	
	// Erase duplicates in arcs
	G.A.assign(se.begin(),se.end());
	G.Narcs = G.A.size();
	cout << "---------------------------------------" << endl;
	
    // return the Graph object
    return G;
}

mycallbackT1::mycallbackT1(const Graph& xG, const vector<GRBVar>& xfa, const vector<int>& xdemand, const int& xUB) {
    // initialize the Callback object
    this->Ncuts = 0;  // the number of cuts, initially 0
    this->G = xG;     // the arcflow graph
    this->fa = xfa;   // the flow
	this->demand = xdemand;
	this->UB = xUB;
}

int mycallbackT1::getNcuts() {
    return this->Ncuts; // return the number of cuts
}

void mycallbackT1::callback() {
    if (where == GRB_CB_MIPSOL) { // if a new incumbent is found
        // find the incumbent solution
        vector<int> f(G.Narcs);                			
        for (int i = 0; i < G.Narcs; i++) f[i] = ceil(getSolution(fa[i]) - EPSILON); 
		
        // find all cycles
		vector<vector<int>> arcsUsed (G.Nvert);
		vector<vector<int>> cycles; 
		int nbRem = 0;
		for (int i = 0; i < G.Narcs; i++) { 
			for (int k = 0; k < f[i]; k++){  
				arcsUsed[G.A[i][0]].push_back(i);
				nbRem++;
			}
		}

		while (nbRem > 0){
			vector<int> cycle;
			int tail = 0; 
			while(arcsUsed[tail].size() == 0) tail++;
			// Create cycle
			while(arcsUsed[tail].size() > 0){
				cycle.push_back(arcsUsed[tail].back()); 
				arcsUsed[tail].pop_back(); 
				tail = G.A[cycle.back()][1];
				nbRem--;
			}
			// Extend cycle
			for(int i = 0; i < cycle.size();i++){
				if(arcsUsed[G.A[cycle[i]][0]].size() > 0){
					tail = G.A[cycle[i]][0];
					int tempI = i;
					while(arcsUsed[tail].size() > 0){
						cycle.insert(cycle.begin() + tempI,arcsUsed[tail].back()); tempI++;
						arcsUsed[tail].pop_back(); 
						tail = G.A[cycle[tempI-1]][1];
						nbRem--;
					}
				}
			}
			cycles.push_back(cycle);
		}

		for(int i = 0; i < cycles.size();i++){
			cout << "Cycle " << i << ": "; 
			for(int j = 0; j < cycles[i].size();j++){
				cout << "( " << G.ITN[G.A[cycles[i][j]][0]][0] << "," << G.ITN[G.A[cycles[i][j]][0]][1] << " " << G.ITN[G.A[cycles[i][j]][1]][0] <<  "," << G.ITN[G.A[cycles[i][j]][1]][1] << " " <<  "," << G.A[cycles[i][j]][2] << ") ";
			}
			cout << endl;
		}

        // for each supplementary cycle, add a no-good cut
		for(int i = 1; i < cycles.size();i++){		
			int count = 0;
			GRBLinExpr LHS = 0; GRBLinExpr RHS = 0;
			vector<bool> isNodeInTheCut (G.Nvert,false);			
			for(int j = 0; j < cycles[i].size();j++){
				isNodeInTheCut[G.A[cycles[i][j]][0]] = true; 
			}
			for (int j = 0; j < G.Narcs; j++) { 
				if (isNodeInTheCut[G.A[j][0]] && !isNodeInTheCut[G.A[j][1]])
					RHS += fa[j];
				if(isNodeInTheCut[G.A[j][0]] && isNodeInTheCut[G.A[j][1]]){
					LHS += fa[j]; 
					if (G.A[j][2] != -1) count += demand[G.A[j][2]];
				}
			}

			// add the cut to the set of variables
            addLazy(LHS <= RHS * (count + UB));  		
            Ncuts++;     
		}
		cout << "Added " << Ncuts << " cuts " << endl;
		cout << "-------------------------------------------------------" << endl;
    }
}

Solution EXACTV1(const Instance& inst) {

	// local variable
	double start = getCPUTime(); 			// starting time
	Solution sol;	sol.Ncuts = 0;			// solution to return
	Graph G = makeGraph1(inst);   
	
    // create a model
    GRBEnv env = GRBEnv();              	// create an environment
    GRBModel model = GRBModel(env);         // create a new model

    // declaration of the variables for the model
	vector<GRBVar> fa(G.Narcs);
	
    // initizalization of the variables for the model
    for (int i = 0; i < G.Narcs; i++) { 
	//	cout << G.ITN[G.A[i][0]][0] << "," << G.ITN[G.A[i][0]][1] << " " << G.ITN[G.A[i][1]][0] << "," << G.ITN[G.A[i][1]][1] << " " << G.A[i][2] << endl;
		if(G.A[i][2] != -1) fa[i] = model.addVar(0, inst.BCs[G.A[i][2]][2], 0, GRB_INTEGER);
		else fa[i] = model.addVar(0, inst.UB, 0, GRB_INTEGER);
    }
    model.update();

	// declare linear expressions
    vector<GRBLinExpr> fIn (G.Nvert, 0);    		// the amount of flow entering each vertex
    vector<GRBLinExpr> fOut (G.Nvert, 0);   		// the amount of flow leaving each vertex
    vector<GRBLinExpr> typeUsed(inst.m + 1, 0); 	// the amount of arcs used of each BC type 
    
	// calculate the linear expressions	
    for (int i = 0; i < G.Narcs; i++) {   		// loop over all arcs
        fIn[G.A[i][1]] += fa[i];        		// inflow
        fOut[G.A[i][0]] += fa[i];       		// outflow
        typeUsed[G.A[i][2]+1] += fa[i];      	// number of BCs used of certain type
    }
	model.update();

    // create flow conservation constraints 
    for (int v = 0; v < G.Nvert; v++) {       	// loop over all vertices
        model.addConstr(fIn[v] == fOut[v]); 	// inflow = outflow
    }

    // create constraint ensuring at least one cycle containing (0, 0)
    model.addConstr(fOut[0] >= 1); 

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
	vector<int> demands; for(int i = 0; i < inst.m; i ++) demands.push_back(inst.BCs[i][2]);
    mycallbackT1 cb = mycallbackT1(G, fa, demands, inst.UB);
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
		
        // find all cycles
		vector<vector<int>> arcsUsed (G.Nvert);
		vector<vector<int>> cycles; 
		int nbRem = 0;
		for (int i = 0; i < G.Narcs; i++) { 
			for (int k = 0; k < ceil(fa[i].get(GRB_DoubleAttr_X) - EPSILON); k++){  
				arcsUsed[G.A[i][0]].push_back(i);
				nbRem++;
			}
		}

		while (nbRem > 0){
			vector<int> cycle;
			int tail = 0; 
			while(arcsUsed[tail].size() == 0) tail++;
			// Create cycle
			while(arcsUsed[tail].size() > 0){
				cycle.push_back(arcsUsed[tail].back()); 
				arcsUsed[tail].pop_back(); 
				tail = G.A[cycle.back()][1];
				nbRem--;
			}
			// Extend cycle
			for(int i = 0; i < cycle.size();i++){
				if(arcsUsed[G.A[cycle[i]][0]].size() > 0){
					tail = G.A[cycle[i]][0];
					int tempI = i;
					while(arcsUsed[tail].size() > 0){
						cycle.insert(cycle.begin() + tempI,arcsUsed[tail].back()); tempI++;
						arcsUsed[tail].pop_back(); 
						tail = G.A[cycle[tempI-1]][1];
						nbRem--;
					}
				}
			}
			cycles.push_back(cycle);
		}

		for(int i = 0; i < cycles.size();i++){
			cout << "Cycle " << i << ": "; 
			for(int j = 0; j < cycles[i].size();j++){
				cout << "( " << G.ITN[G.A[cycles[i][j]][0]][0] << "," << G.ITN[G.A[cycles[i][j]][0]][1] << " " << G.ITN[G.A[cycles[i][j]][1]][0] << "," << G.ITN[G.A[cycles[i][j]][1]][1] << " " <<  "," << G.A[cycles[i][j]][2] << ") ";
			}
			cout << endl;
		}
		
		if(cycles.size() == 1){
			cout << "No subtour, solution is optimal" << endl;
			// get bin for each BC type
			sol.assignedBin.resize(inst.m);
			int binIdx = 0;
			for(int i = 0; i < cycles[0].size();i++){
				if(G.A[cycles[0][i]][2] >= 0) sol.assignedBin[G.A[cycles[0][i]][2]].push_back(binIdx);
				else binIdx++;
			}
			if(sol.LB == sol.UB) sol.opt = 1;
		}
    }

    sol.timeT = getCPUTime() - start;	
    return sol;
}
