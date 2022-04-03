#include "main.h"

int main(int argc, char **argv){          

	// Read input and output paths
	string path = argv[1];	
	string filein = argv[2];
	string pathAndFileout = argv[3];
	
    // initialize the input variables from a file
    Instance inst = readInstance(path + filein);
	inst.print();
   
    // find the NAIVE solution
    Solution sol_NAIVE = NAIVE(inst);
	printInfo(pathAndFileout, sol_NAIVE, filein);
	printSInfo(sol_NAIVE, inst);
	
}

