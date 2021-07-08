#ifndef Pop_H
#define Pop_H

#include <iostream>
#include <string>
#include <stdlib.h>
#include "Bac.h"
#include <vector>
#include <memory>
#include <unordered_map>

using namespace std;

//This class represents one population and handles most of the operations of this population. 

class Pop {
	
public:
	int popindex; // Which population this one is in the Superpop

	vector<shared_ptr<Bac>> bacs; // Vector to store the Bac objects
	//shared_ptr's are used for the superior efficiency and usability.

	Para parand; // Para object for this Pop object to use. 

	int ngenes;
	int nbacs;
	int genelength;
	int tmuts; 
	int snpsfromuts; 
	int mutsonsnp; 
	int mutstoanc; 
	int nrecombinations;
	int nrectries;
	int nullrecombs;
	int recsnps; 
	int addrecn;  
	int dsnps;
	int delrecn;
	double mutlambda;
	double reclambda;
	double inslambda;
	double dellambda;

	map<int,double> mismadistr;

	double GCs; // Currently not in use
	int aa; // Amounts of bases. Currently not in use
	int cc;		
	int gg;	
	int tt;
	int bacsimmigrated;
	int bacsemigrated;
	int nmepidems;
	int nmpedimbacs;
	int miglambda;
	string  mismastr;
	string  allemismastr;
	string  Popmismastr;
	string  summarystr;

	Pop(int nrbacs, int nrgenes, int genelen, Para pa);
	// Constructs population with clones of new bacterium with randomly generated sequences as loci.


	Pop(int popi, int nbacteria, shared_ptr<Bac> bacinput, Para pa); 
	// Constructs population based bacterium object (Bac) given as argument.
	
	void printPop(int interval);
	void mutate(int nmuts, int thisgeneration); 
	int recombineSegment(int nrecs) ; // Recombinations on random segments of loci
	int recombineGene(int nrecs) ; // Recombinations on random whole loci 
	int insertions(int nins) ;
	int deletions(int ndel) ; 
	string migrate(Pop* dest, int migamount); // Prompts migration from this population to dest.
	string mepidemic(); 
	void select(); 
	void GeneWiseDiffs(const vector<int>  & randbacinds); 
	void InterpopGeneWiseDiffs(const vector<int>  & randbacinds,const vector<int>  & Brandbacinds,const Pop & popB); 
	double Burner(const vector<int>  & randbacinds); 
	void GeneDistance(const shared_ptr<Gene> & igene, const shared_ptr<Gene> & jgene, int & alldifs, int & sharedlen, int & mutationlydifs); 
	// Compares the given loci of two bacs
	string summary(int generation, int popindex); // Prints basic summary information of population.
	void PopInfo();
	void BacRangeInfo(int lowbound, int upbound); 
	// Prints some of the data about the bacteria in the given range from population.
	void BacRangeInfo(int lowbound, int upbound, int m); 
	// Prints some of the data about the specified gene m of the bacteria in the given range from population.
	int Diversity(int generation, int popindex);

};


#endif
