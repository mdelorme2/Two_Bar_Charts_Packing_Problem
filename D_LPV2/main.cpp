#include "main.h"

int main(int argc, char **argv){          

	// Read input and output paths
	string path = argv[1];	
	string filein = argv[2];
	string pathAndFileout = argv[3];
	
    // initialize the input variables from a file
    Instance inst = readInstance(path + filein); 
	inst.print();
   
    // find the LPV2 solution
    Solution sol_LPV2 = LPV2(inst);
    printInfo(pathAndFileout, sol_LPV2, filein);
}

