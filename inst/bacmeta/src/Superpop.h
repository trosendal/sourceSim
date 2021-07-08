#ifndef Superpop_H
#define Superpop_H
#include "Para.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include "Pop.h"
#include <vector>
#include <unordered_map>


using namespace std;
//This class stores the populations and handles most of the top, i.e. metapopulation, level operations and prompts lower level operations.
class Superpop {
	
	
public:
	Para* parand; // Used for RNGs and parameters in Superpop functions.
	string fnadd;
	int npops;
	
	vector<Pop> Pops; // Stores Population objects
	
	static int thisgeneration; // Stores current generation

	int nrecs; 
	int accrecs; 
	int ntotalmuts; 
	int ntotalaccrecs; 
	int nmigevent; 
	
	vector<vector<double> > migmatrix; /* Vector of vectors, i.e. matrix, for migration parameters between pop_i and pop_j. */
		
	string mics;
	string migs;
	bool miginit;

	bool distinit;
	
	vector<vector<int> > distmatrix; /* Vector of vectors, i.e. matrix, of sample sizes for interpopulation distance computing between pop_i and pop_j. */
	int maxdistsample;

	static int allcount;
	static int stcount;
	static int backtoinitialmuts;
	
	static vector<unordered_map<string, int>> SeqTypes; 
	static unordered_map<string, int> MLSTs;
	
	static map<int,int> recs; //<HD/SimilarityTestFail,#SuchEvents> Similarity test fail as -1
	
	static map<int,int> indels;
		
	Superpop(int nrpops, int nrbacs, int nrgenes, int genelen, Para pa); 
	~Superpop(); 

	void printSuperpop(int bacinterval, int popinterval); 
	
	void startmutations();
	
	void startrecombinations();
	
	void startindels();
	
	bool Migmatrixupd();
	
	void migration();
	
	void mepidemic();
	
	void selection();
	
	void newGeneration();
	//Runs mutation, indels recombination, migration and selection methods. Updates generation count
	
	bool distmatrixupd();

	static int getAllele(string str, int gene);
	static int getMLST(vector<shared_ptr<Gene>> &genes);
	
	void document();
	void savedocumentation();
	


};

#endif
