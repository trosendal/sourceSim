#ifndef Bac_H
#define Bac_H

#include <iostream>
#include <string>
#include <stdlib.h>
#include "Gene.h"
#include <vector>
#include <map>
#include <memory>
#include <set>	


using namespace std;

//This class represents one bacterium and runs the events of the bacterium. 

class Bac {
	
public:

	vector<shared_ptr<Gene>> genes; // Vector of shared_ptrs to Gene objects that represent the genes.

	int ngenes; 
	int snpn;
	int emergeneration;
	vector<int> mepis;
	int genelength; 
	int nmutations; 
	int nrecombinations; 
	int MLST;
	int mlstdo; 
	int recsitesdo;

	double GC; // Variable for keeping track of GC ratio. Currently not in use
	int aa; // Amounts of bases. Currently not in use
	int cc;		
	int gg;	
	int tt;

	
	
	Bac(int ngenes, int genelen, Para pa); 
	
	void printGenes(); 
	void mutate(Para &pa, int& snpd); 
	void recombine(int recgene, int rocus, string donatedseq, int snpchange, Para& pa); 
	void insertion(Para& pa);
	void deletion(Para& pa);
	vector<double> getbases(int curgeneration, Para& pa); // Currently not in use

};


#endif
