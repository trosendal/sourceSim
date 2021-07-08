#ifndef Gene_H
#define Gene_H
#include "Para.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>

using namespace std;
class Gene {
	
public:
	
	string genstr;	// String of the main sequence, i.e. not including insertions

	int genei; //  Which loci this one is.


	double GC; //Currently not in use.
	int aa; // Variables for amounts of bases. Currently not in use.
	int cc;		
	int gg;	
	int tt;
	int genelength;
	int nmutations;
	int nrecombinations;

	int updategeneration; 

	int allele;
	int mlstdo;
	int recsitesdo;
	map<int, char> snps;
	
	map<int, string> ins;
	map<int, vector<pair<int,int>>> insrecs;
	int inslens;
	
	set<int> dels;
	map<int, int> recsites; //list of <endpoint,startpoint> -pairs representing recombination sites. considering only recombinations with changes in snp count


	static double mutmatrix[12];

	Gene(int length, Para pa, int ng); 
	string getGeneText(bool inserts);  
	void mutate(Para &pa, int& snpd);
	char mutationBase(char oldbase, double baserand);
	void recombine(int rocus, string donatedseq); 
	void addRecsite (int endpoint, int startpoint);
	void mutationRecsitesplit(int position);	
	void mutationInsRecsitesplit(int inskey, int position);	
	int recsiteSNPs(int inskey); 
	void insertion(Para &pa);	
	void deletion(int rocus, int len);
	void updateBasecounts();	
	void updbasecountbySNP(int snpposition, char snpbase);
};


#endif
