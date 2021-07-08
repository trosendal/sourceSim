#ifndef PARA_H
#define PARA_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <stdio.h>
#include <random>


using namespace std;
// Header file of Para class. Used to declare variables for the Para object.

struct foo // Structure object to store RNGs
{
	static mt19937 generator;
	static mt19937 generator2;
	static uniform_real_distribution<double> unidistr;
	static poisson_distribution<int> poisdistr;
	static geometric_distribution<int> geodistr;
	
};


class Para {
public:
	
	static string rands;
	static string drands;
	static int ndrands;
	static int staticseed;
	static int trueseed;
	
	static string inputfilename;
	static string migfilename;
	static string outmod;
	//For zipf
	static long double qd,iqd,qi,iqi,dH0,ds,dHmax,iH0,is,iHmax;
	

	// Input parameters for the Para object. Add new parameters here according to previous ones	
	string fna;  // Filename modifier 
	int seed; // Set seed or give zero to have random seed from time. 
	int ngenerations; //number of generations to run
	double domean; //Compute and save mean distance summary every 1000 generation (1) or not (0)
	int genelen; //length of a gene
	int ngenes; //number of genes in bacteria
	int nbacs; //number of bacteria in population
	int npops; //number of populations in Superpop
	int debugval; //interval of debug
	int debugpop; //debug population
	int custommigs; // switch for migration matrix input
	double migscaler; //migration rate scaler
	double migProp; //migration probability
	double mepiAmo; //microepidemic probability
	double mepiSize; //mean microepidemic size scaler
	double mutrate; //mutation rate
	double insrate; //insertion rate
	double insa; //insertion length parameter
	double delrate; //deletion rate
	double dela; //deletion length parameter
	double indelmax; //deletion length parameter
	double recrate; //recombination rate
	double reclen; //recombination mean length
	double recacc; //recombination acceptance parameter
	double recsites; //recombination sites gathering
	double aProp,tProp,gProp,cProp; //proportion of bases in gene initialization
	double atot,atog,atoc; //probability of mutation from base x to y
	double ttoa,ttog,ttoc;
	double gtoa,gtot,gtoc;
	double ctoa,ctot,ctog;
	int sele; // keep population as is (0) or randomly select for new generation (1)
	int rande; // keep evolution event ordered the same (0) or randomize each generation (1)
	int summaryinterval; // interval of generations for summary printing
	double genewisemisma; // Count lociwise mismatches with this proportion
	double mutonly; // Count mutation-only lociwise mismatches (1) or not (0)
	int genewiseval;	//	Document gene by lociwise mismatces with this generation interval
	double popmisma; // Count interpopulation lociwise mismatches with this proportion
	double popmutonly; // Count interpopulation mutation-only lociwise mismatches  (1) or not (0)
	int popeval;	//	Document interpopulation lociwise mismatces with this generation interval
	double iniseq; // Save initial genome (1) or not (0)
 	double seqsample; // Proportion of population for saving sequences
	int seqval; // Save sequences with this generation interval
	int diversityinterval; // Save diversity with this interval
	int mutsinterval; //Do mutation documenting (1) or not (0)
	int recsinterval; //Do recombination documenting with this generation interval
	int recdocthresh; //Do recombination event documentation with this hd threshold
	int indeldo; //Do indel documenting (1) or not (0)
	int migdo; //Do migration documenting (1) or not (0)
	int micdo; //Do microepidemic documenting (1) or not (0)
	int mlstdo; //Do MLST documenting (1) or not (0)
	int alleseq; //Do sequence-allele type documenting (1) or not (0)
	
	
	//Variables for Para object
	bool selerand; // boolean to store input from sele
	bool randevents; // boolean to store input from randevents
	vector<double> baseProps;
	static vector<string> origenes;
	foo ranstruct;
	
	double geop; // Probability for recombination lenght geometric distribution
	
	//For zipf
	long double U,X,K;
	#define H(x,q,iq,v) exp(q*log(x+v))*iq
	#define H1(x,q,iq,v) -v+exp(iq*log(q*x))	
	
	//Functions of Para object
	Para(); // Constructor of Para for creating Para object
	~Para(); // Destructor of Para object
	vector<double> getBaseprops(); // Gives the distribution of bases in a vector 
	int getRandom(); //To get a random integer from uniform distribution [0, 2^32-1]
	double getZeroone(); //To get a random number from uniform distribution [0,1]
	double getUtilZeroone(); //To get a random number from uniform distribution [0,1]
	vector<double> getZOvector(int amount); //To get a random number from uniform distribution [0,1]
	vector<int> getRandvector(int amount, int upperbound); //To get a random number from uniform distribution [0,1]
	vector<int> getUtilNoRepSample(const vector<int> & presample, int amount); // To get random sample without replacement for documentation purposes. The sample is either from  range [0, presample[0]-1] if size of presample vector is 1, or otherwise sample of the values in presample.
	int getPoisrand(double mean); // To get a random number from Poisson distribution with expected value of mean
	int getGeorand(); // To get a random number from Poisson distribution with expected value of mean
	void initZipf(double vv); //Set up Zipf distribution random number generation, as in  Hörmann, W. and G. Derflinger, G. (1996)
	int getZipf(bool iord, double vv); // To get a random number from Zipf distribution with parameter a. iord true for deletion, false for insertion.	
	void updateSeed();  // To update the seed of the random generator
	int readInput(bool defaultfile); // Input file reader 
	static int parser(char inputfn[80] , string outputfn); // Converts old style parameter input file to new format based on parameter codes. Only here for legacy reasons.

};
#endif
