#include "Para.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <stdio.h>
#include <random>
#include <chrono>
#include <limits>

using namespace std;
/* Para objects are used to read parameters from simu*.input file, provide them to other classes and generate random numbers from different distributions. 
*/


uniform_real_distribution<double> foo::unidistr(0.0, 1.0);  
mt19937 foo::generator(0);  // RNGs for simulation core
mt19937 foo::generator2(0);  // RNGs for documentation sampling
string Para::rands="";
string Para::drands="";
int Para::ndrands=0;
int Para::staticseed=0;
int Para::trueseed=0;
long double Para::qd=0;  // These long doubles are used for Zipf-distribution
long double Para::iqd=0;
long double Para::qi=0;
long double Para::iqi=0;
long double Para::dH0=0;
long double Para::ds=0;
long double Para::dHmax=0;
long double Para::iH0=0;
long double Para::is=0;
long double Para::iHmax=0;

vector<string> Para::origenes;


Para::Para(){
	
	cout<< "Reading parameter input and initializing randomness." << endl;


	if (readInput(1)==0) { // Read default file, abort if non-successful return.
		outmod="error";
		return; 
	}  
	
	if (readInput(0)==0) { // Read parameter input file, abort if non-successful return.
		outmod="error";
		return;
	}
	
	// fscanf function doesn't support reading boolean so here booleans are initialized. 
	selerand = sele; 
	randevents = rande;
	
	geop=1.0/reclen;
	
	
	if(outmod!="") {   // If output modifier is given with command line flag
		fna=outmod;
	}
	
	if(Para::migfilename!="migration.input") {   // If custom migration input file is given with command line flag, override universal migration rate.
		custommigs=1;
	}
	
	if (mutonly==1) {  // If mutation-only genewise distances are requested, override recsite documenting setting to "on". 
		recsites=1;
	}
	
	
	if ((seed==0) && (Para::trueseed==0) ) {  // If requested, random number generator is seeded with system clock time instead of given seed.
		seed=chrono::system_clock::now().time_since_epoch().count();
		Para::trueseed=seed;		// Static global variables are used to handle multiple instances of Para.
		Para::staticseed=seed;
	}
	if ((seed!=0) && (Para::trueseed==0) ) {
		Para::trueseed=seed;
		Para::staticseed=seed;
	} 
	
	if (mutsinterval==0 ) {
		mutsinterval=recsinterval;
	} 
	
	foo::generator.seed(Para::staticseed); 
	foo::generator2.seed(getZeroone()); 
	Para::staticseed= foo::generator();
	
	baseProps={aProp,aProp+tProp,aProp+tProp+gProp};
}

Para::~Para() { }	



vector<double> Para::getBaseprops() { /* Gives a vector of the proportions of bases for initialization of genes.*/
	return baseProps;
}


int Para::getRandom() {  // Gives a random integer from uniform distribution [0, 2^32-1]
	int randvalue= foo::generator();
	return randvalue;
}

double Para::getZeroone() { // Gives random real number from uniform distribution in [0,1] for simulation purposes
	double randvalue =foo::unidistr(foo::generator);
	return  randvalue;
}

double Para::getUtilZeroone() { // Gives random real number from uniform distribution in [0,1] for documentation purposes
	double randvalue =foo::unidistr(foo::generator2);
	return  randvalue;
}

vector<double> Para::getZOvector(int amount) { // Gives requested amount of random real numbers from uniform distribution in [0,1]

	vector<double> randvector; 
	randvector.reserve(amount);
	for (int i =0; i< amount ; i++) {
		Para::ndrands++;
		randvector.push_back(foo::unidistr(foo::generator));
	}	
	return  randvector;
}

vector<int> Para::getRandvector(int amount, int upperbound) { // Gives requested amount of random real numbers from uniform distribution in [0,upperbound]

	vector<int> randvector; 
	randvector.reserve(amount);
	for (int i =0; i< amount ; i++) {
		randvector.push_back(foo::unidistr(foo::generator)*upperbound);
	}
	return  randvector;
}

vector<int> Para::getUtilNoRepSample(const vector<int> & presample, int amount) { // To get random sample without replacement for documentation purposes. The sample is either from  range [0, presample[0]-1] if size of presample vector is 1, or otherwise sample of the values in presample. Modified from algorithm presented in C. T. Fan, M. E. Muller, and I. Rezucha (1962) and T. G. Jones (1962).

	vector<int> randinds(amount);

    int N;

    if ((int)presample.size() ==1 ) {
    	N=presample[0];
    }
    else {
    	N=(int)presample.size();
    }

    if( N==amount) {
		iota (begin(randinds), end(randinds), 0);
		return  randinds;
	}
	else if (N<amount) {
		cerr << "Invalid sample size: " << amount << ".";
		return vector<int> {};
	}

	int t = 0; 
	int m = 0; 
	double u;

	while (m < amount)
	{
		u = foo::unidistr(foo::generator2); 

		if ( (N - t)*u >= amount - m )
		{
			t++;
		}
		else
		{
   			 if ((int)presample.size() ==1 ) {
		    	randinds[m] = t;
		    }
		    else {
		    	randinds[m] = presample[t];
		    }

			t++;
			m++;
		}
	}	


	return randinds;

}

int Para::getPoisrand(double mean) { // Gives random number from a Poisson distribution with the given mean 
	poisson_distribution<int> poissodistr(mean);
	return poissodistr(foo::generator);
}

int Para::getGeorand() { // Gives random number from Poisson distribution with expected value of 4
	geometric_distribution<int> geodistr(geop);
	return geodistr(foo::generator);
}

void Para::initZipf(double vv) { //Sets up Zipf distribution random number generation, as in  Hörmann, W. and G. Derflinger, G. (1996)
		
		unsigned int imax = ~0;
		double v=vv;
		//For deletions
			//dela=q
			Para::qd=1.0-dela;	
			Para::iqd=1.0/Para::qd;
			Para::dH0=H(0.5,Para::qd,Para::iqd,v)-exp(log(v)*-dela);
			double tempX=H(1.5,Para::qd,Para::iqd,v)-exp(log(v+1)*(-dela));
			Para::ds= 1.0-H1(tempX,Para::qd,Para::iqd,v);
			//cout << "ds " << ds << endl;
			Para::dHmax = H(imax + 0.5 ,Para::qd,Para::iqd,v);
		
		//For insertions			
			Para::qi=1.0-insa;	
			Para::iqi=1.0/Para::qi;
			
			Para::iH0=H(0.5,Para::qi,Para::iqi,v) - exp(log(v)*-insa);
			tempX=H(1.5,Para::qi,Para::iqi,v)-exp(log(v+1.0)*(-insa));
			Para::is= 1.0-H1(tempX,Para::qi,Para::iqi,v);
			//cout << "is " << is << endl;
			Para::iHmax = H(imax + 0.5 ,Para::qi,Para::iqi,v);
			
}

int Para::getZipf(bool iord, double vv) { // Gives random number from Zipf distribution with parameter "a". "iord" is TRUE for deletion, FALSE for insertion.
	double v=vv;
	int limit= 1e6;
	if (iord) { // Zipf for deletion
		   while (limit--) {
			U=getZeroone();
			U=Para::dHmax+(U*(Para::dH0-Para::dHmax));
			X=H1(U,Para::qd,Para::iqd,v);
			K=floor(X+0.5);
			if (K-X<=Para::ds) return int(K+1);
			else if ( U >=H(K+0.5,Para::qd,Para::iqd,v)-exp(-1*log(v+K)*(dela+1)) ) return int(K+1);
		}    
		
		
		
	}	
	else { // Zipf for insertion
		while (limit--) {
			U=getZeroone();
			U=Para::iHmax+(U*(Para::iH0-Para::iHmax));
			X=H1(U,Para::qi,Para::iqi,v);
			K=floor(X+0.5);
			
			if (K-X<=Para::is) return int(K+1);
			else if ( U >=H(K+0.5,Para::qi,Para::iqi,v)-exp(-1*log(v+K)*(insa+1))) return int(K+1);
		}	
	}
	cerr << "Something went wrong with Zipf RNG" << endl;
	return 1;
}


void Para::updateSeed(){ // Updates the seed of the random generator with random number from the generator.
	int seedupdater;
	seedupdater=foo::generator();
	foo::generator.seed(seedupdater);
}


int Para::readInput(bool defaultfile) {  // Reads parameter inputs based on parameter codes.



	char filename[80];
	char dummy[180];
	char tempstringreader[80];
	

	if (defaultfile) { 
		strcpy(filename,"default.params");
	}
	else {
		strcpy(filename,Para::inputfilename.c_str());
	}
	
	
	FILE* infile;
	infile = fopen(filename,"r");

	char paramcode[5];
	strcpy(paramcode,"null");

	if (infile != NULL) {	

		int breaker=0;

		while (strcmp(paramcode,"#")!=0 && breaker<123) {

			breaker++;

			fscanf(infile,"%s",paramcode);

			if(strcmp(paramcode,"OPFN:")==0) {
				fscanf(infile,"%s",tempstringreader);fgets (dummy, 180, infile);
				fna=string(tempstringreader);
				}
			else if (strcmp(paramcode, "SEED:")==0) { 
				fscanf(infile,"%d",&seed);fgets (dummy, 180, infile);
				} 
			else if (strcmp(paramcode, "GENR:")==0) { 
				fscanf(infile,"%d",&ngenerations);fgets (dummy, 180, infile);		
				}
			else if (strcmp(paramcode, "MEAN:")==0) { 
				fscanf(infile,"%lf",&domean);fgets (dummy, 180, infile);
				}
			else if (strcmp(paramcode, "LOLE:")==0) { 
				fscanf(infile,"%d",&genelen);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "NLOC:")==0) { 
				fscanf(infile,"%d",&ngenes);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "NBAC:")==0) { 
				fscanf(infile,"%d",&nbacs);fgets (dummy, 180, infile);						
				
				}
			else if (strcmp(paramcode, "NPOP:")==0) { 
				fscanf(infile,"%d",&npops);fgets (dummy, 180, infile);			
				
				}
			else if (strcmp(paramcode, "PRIG:")==0) { 
				fscanf(infile,"%d",&debugval);fgets (dummy, 180, infile);			
				
				}
			else if (strcmp(paramcode, "PRIP:")==0) { 
				fscanf(infile,"%d",&debugpop);fgets (dummy, 180, infile);	
				
				}
			else if (strcmp(paramcode, "MIGI:")==0) { 
				fscanf(infile,"%d",&custommigs);fgets (dummy, 180, infile);		
				
				}
			else if (strcmp(paramcode, "MIGR:")==0) { 
				fscanf(infile,"%lf",&migscaler);fgets (dummy, 180, infile);		
				}
			else if (strcmp(paramcode, "MIGP:")==0) { 
				fscanf(infile,"%lf",&migProp);fgets (dummy, 180, infile);		
				
				}
			else if (strcmp(paramcode, "MICA:")==0) { 
				fscanf(infile,"%lf",&mepiAmo);fgets (dummy, 180, infile);		
				
				}
			else if (strcmp(paramcode, "MICS:")==0) { 
				fscanf(infile,"%lf",&mepiSize);fgets (dummy, 180, infile);		
				
				}
			else if (strcmp(paramcode, "MUTR:")==0) { 
				fscanf(infile,"%lf",&mutrate);fgets (dummy, 180, infile);		
				
				}
			else if (strcmp(paramcode, "RECR:")==0) { 
				fscanf(infile,"%lf",&recrate);fgets (dummy, 180, infile);		
				
				}
			else if (strcmp(paramcode, "RECL:")==0) { 
				fscanf(infile,"%lf",&reclen);fgets (dummy, 180, infile);		
				
				}
			else if (strcmp(paramcode, "RECA:")==0) { 
				fscanf(infile,"%lf",&recacc);fgets (dummy, 180, infile);		
				
				}
			else if (strcmp(paramcode, "RECS:")==0) { 
				fscanf(infile,"%lf",&recsites);fgets (dummy, 180, infile);	
				
				}
			else if (strcmp(paramcode, "INSR:")==0) { 
				fscanf(infile,"%lf",&insrate);fgets (dummy, 80, infile);		
				
				}
			else if (strcmp(paramcode, "INSL:")==0) { 
				fscanf(infile,"%lf",&insa);fgets (dummy, 80, infile);		
				
				}
			else if (strcmp(paramcode, "DELR:")==0) { 
				fscanf(infile,"%lf",&delrate);fgets (dummy, 80, infile);		
				
				}
			else if (strcmp(paramcode, "DELL:")==0) { 
				fscanf(infile,"%lf",&dela);fgets (dummy, 80, infile);		
				
				}
			else if (strcmp(paramcode, "INDM:")==0) { 
				fscanf(infile,"%lf",&indelmax);fgets (dummy, 80, infile);		
				
				}
			else if (strcmp(paramcode, "PROA:")==0) { 
				fscanf(infile,"%lf",&aProp);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "PROT:")==0) { 
				fscanf(infile,"%lf",&tProp);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "PROG:")==0) { 
				fscanf(infile,"%lf",&gProp);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "PROC:")==0) { 
				fscanf(infile,"%lf",&cProp);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "ATOT:")==0) { 
				fscanf(infile,"%lf",&atot);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "ATOG:")==0) { 
				fscanf(infile,"%lf",&atog);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "ATOC:")==0) { 
				fscanf(infile,"%lf",&atoc);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "TTOA:")==0) { 
				fscanf(infile,"%lf",&ttoa);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "TTOG:")==0) { 
				fscanf(infile,"%lf",&ttog);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "TTOC:")==0) { 
				fscanf(infile,"%lf",&ttoc);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "GTOA:")==0) { 
				fscanf(infile,"%lf",&gtoa);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "GTOT:")==0) { 
				fscanf(infile,"%lf",&gtot);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "GTOC:")==0) { 
				fscanf(infile,"%lf",&gtoc);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "CTOA:")==0) { 
				fscanf(infile,"%lf",&ctoa);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "CTOC:")==0) { 
				fscanf(infile,"%lf",&ctot);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "CTOG:")==0) { 
				fscanf(infile,"%lf",&ctog);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "RSEL:")==0) { 
				fscanf(infile,"%d",&sele);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "RORD:")==0) { 
				fscanf(infile,"%d",&rande);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "SUMI:")==0) { 
				fscanf(infile,"%d",&summaryinterval);fgets (dummy, 180, infile);
				
				}
			else if (strcmp(paramcode, "GWDS:")==0) { 
				fscanf(infile,"%lf",&genewisemisma);fgets (dummy, 180, infile);
 				
				}
			else if (strcmp(paramcode, "MGWD:")==0) { 
				fscanf(infile,"%lf",&mutonly);fgets (dummy, 180, infile);
 			
				}
			else if (strcmp(paramcode, "GDWI:")==0) { 
				fscanf(infile,"%d",&genewiseval);fgets (dummy, 180, infile);
 
				}
			else if (strcmp(paramcode, "GWDI:")==0) { 
				fscanf(infile,"%d",&genewiseval);fgets (dummy, 180, infile);
 
				}
			else if (strcmp(paramcode, "PGWS:")==0) { 
				fscanf(infile,"%lf",&popmisma);fgets (dummy, 180, infile);
 				
				}
			else if (strcmp(paramcode, "PMGW:")==0) { 
				fscanf(infile,"%lf",&popmutonly);fgets (dummy, 180, infile);
 			
				}
			else if (strcmp(paramcode, "PGWI:")==0) { 
				fscanf(infile,"%d",&popeval);fgets (dummy, 180, infile);
 			
				}
			else if (strcmp(paramcode, "ISEQ:")==0) { 
				fscanf(infile,"%lf",&iniseq);fgets (dummy, 180, infile);
 
				}
			else if (strcmp(paramcode, "SEQS:")==0) { 
				fscanf(infile,"%lf",&seqsample);fgets (dummy, 180, infile);
 
				}
			else if (strcmp(paramcode, "SEQI:")==0) { 
				fscanf(infile,"%d",&seqval);fgets (dummy, 180, infile);
 
				}
			else if (strcmp(paramcode, "STRA:")==0) { 
				fscanf(infile,"%d",&diversityinterval);fgets (dummy, 180, infile);
 
				}
			else if (strcmp(paramcode, "MUTD:")==0) { 
				fscanf(infile,"%d",&mutsinterval);fgets (dummy, 180, infile);
 
				}
			else if (strcmp(paramcode, "RECI:")==0) { 
				fscanf(infile,"%d",&recsinterval);fgets (dummy, 180, infile);
 
				}
			else if (strcmp(paramcode, "RECT:")==0) { 
				fscanf(infile,"%d",&recdocthresh);fgets (dummy, 180, infile);
 
				}
			else if (strcmp(paramcode, "INDD:")==0) { 
				fscanf(infile,"%d",&indeldo);fgets (dummy, 80, infile);
 
				}
			else if (strcmp(paramcode, "MIGD:")==0) { 
				fscanf(infile,"%d",&migdo);fgets (dummy, 180, infile);
 
				}
			else if (strcmp(paramcode, "MICD:")==0) { 
				fscanf(infile,"%d",&micdo);fgets (dummy, 180, infile);
 
				}
			else if (strcmp(paramcode, "MLST:")==0) { 
				fscanf(infile,"%d",&mlstdo);fgets (dummy, 180, infile);
 
				}
			else if (strcmp(paramcode, "STAL:")==0) { 
				fscanf(infile,"%d",&alleseq);fgets (dummy, 180, infile);
 
				}
			else if (strcmp(paramcode, "#")==0) { 
				cout << "Parameter reading completed. " << endl;
			}
			
			else {
				cerr << "Unknown parameter code. Aborting." << endl;
				return 0;
			}
			 
		}	


	} else {
		cerr << "Error while opening input file: '" << filename << "'. Aborting. " << endl; 

		return 0;
	}
	fclose(infile);



	if (defaultfile) { 
		cout << "Parameters initialized." << endl; 
	}
	else {
		cout << "User input parameters set based on file: '"<< filename << "'." << endl; 
	}

	

		

	return 1;

}

 
int Para::parser(char inputfn[80] , string outputfn) { // Converts old style parameter input file to new format based on parameter codes. Only here for legacy reasons.


	FILE* infile;
	infile = fopen(inputfn,"r");


	char value[80];
	char variable[180];
	strcpy(value,"null");


	string parsed="";

	cout << "Beginning converting: \n" << endl;

	if (infile != NULL) {

		int breaker=0;

		while (strcmp(value,"#")!=0 && breaker<123) {

			breaker++;

			fscanf(infile,"\n%s",value);
			// fscanf(infile,"%s",variable);
			fgets(variable,160,infile);

			// cout << value << "\t\t" << variable << " |" << endl;

			char* temp=variable;

    		while(*temp && isspace(*temp)) temp++;



			cout << value << "\t\t" << temp ;

			if(strncmp(temp, "#Output file name modifier",20) ==0) {
				parsed+="OPFN: ";
				parsed+=value;
				parsed+="\t\t";
				parsed+=temp;
			}
			else if(strncmp(temp, "#Set seed or give zero",20) ==0) {
				parsed+="SEED: ";
				parsed+=value;
				parsed+="\t\t";
				parsed+=temp;
			}
			else if(strncmp(temp, "#Number of generations to run",20) ==0) {
				parsed+="GENR: ";
				parsed+=value;
				parsed+="\t\t";
				parsed+=temp;
			} 
			else if(strncmp(temp, "#Compute and save mean distance",20) ==0) {
				parsed+="MEAN: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
			}
			else if(strncmp(temp, "#Deviance limit for convergence acceptance :=",20) ==0) {
				parsed+="MEAN: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
			}
			else if(strncmp(temp, "#Length of a gene",15) ==0) {
				parsed+="LOLE: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
			}
			else if(strncmp(temp, "#Number of genes in bacteria",15) ==0) {
				parsed+="NLOC: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
			}
			else if(strncmp(temp, "#Number of bacteria in",15) ==0) {
				parsed+="NBAC: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
			}
			else if(strncmp(temp, "#Number of populations in",15) ==0) {
				parsed+="NPOP: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Debugging output generation, ",25) ==0) {
				parsed+="PRIG: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Debugging output population,",25) ==0) {
				parsed+="PRIP: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Using (1) migration rates",15) ==0) {
				parsed+="MIGI: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Mean migration rate scaler",15) ==0) {
				parsed+="MIGR: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Migration probability",15) ==0) {
				parsed+="MIGP: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Mean microepidemic amount scaler",25) ==0) {
				parsed+="MICA: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Mean microepidemic size scaler",25) ==0) {
				parsed+="MICS: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Mutation rate per nucleotide",20) ==0) {
				parsed+="MUTR: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Recombination rate in relation ",30) ==0) {
				parsed+="RECR: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Recombination mean length. Give ",30) ==0) {
				parsed+="RECL: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Recombination acceptance parameter ",30) ==0) {
				parsed+="RECA: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Gather recombination site metadata",30) ==0) {
				parsed+="RECS: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#insertion rate in relation to mutations",30) ==0) {
				parsed+="INSR: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#insertion length parameter for Zipf",30) ==0) {
				parsed+="INSL: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#deletion rate in relation to",20) ==0) {
				parsed+="DELR: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#deletion length parameter for Zipf ",20) ==0) {
				parsed+="DELL: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#maximum indel length as a",20) ==0) {
				parsed+="INDM: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Proportion of base A",21) ==0) {
				parsed+="PROA: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Proportion of base T",21) ==0) {
				parsed+="PROT: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Proportion of base G",21) ==0) {
				parsed+="PROG: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Proportion of base C",21) ==0) {
				parsed+="PROC: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#P(a->t) ",8) ==0) {
				parsed+="ATOT: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#P(a->g)",8) ==0) {
				parsed+="ATOG: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#P(a->c)",8) ==0) {
				parsed+="ATOC: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#P(t->a)",8) ==0) {
				parsed+="TTOA: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#P(t->g)",8) ==0) {
				parsed+="TTOG: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#P(t->c)",8) ==0) {
				parsed+="TTOC: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#P(g->a)",8) ==0) {
				parsed+="GTOA: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#P(g->t)",8) ==0) {
				parsed+="GTOT: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#P(g->c)",8) ==0) {
				parsed+="GTOC: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#P(c->a)",8) ==0) {
				parsed+="CTOA: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#P(c->t)",8) ==0) {
				parsed+="CTOC: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#P(c->g)",8) ==0) {
				parsed+="CTOG: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Randomly select bacteria",20) ==0) {
				parsed+="RSEL: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Keep order of evolution ",20) ==0) {
				parsed+="RORD: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Summary info: documentation ",20) ==0) {
				parsed+="SUMI: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Gene-wise mismatches and mutation-only mismatches: Size as proportion of the population.", 87) ==0) {
				parsed+="GWDS: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Do mutation-only gene-wise mismatches",20) ==0) {
				parsed+="MGWD: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Gene-wise mismatches and mutation-only mismatches: Generation interval. Documented on-the-go to generation specific files.", 90) ==0) {
				parsed+="GWDI: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Sequence saving: Sample size as proportion",30) ==0) {
				parsed+="SEQS: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Sequence saving: Generation interval",30) ==0) {
				parsed+="SEQI: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Strain-id composition documenting interval",30) ==0) {
				parsed+="STRA: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Do mutation documenting with this interval,",30) ==0) {
				parsed+="MUTD: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Do recombination documenting with this interval",30) ==0) {
				parsed+="RECI: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Do recombination event documenting with this hd threshold or not (0).",30) ==0) {
				parsed+="RECT: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Do indel documenting pp(1) or not (0)",30) ==0) {
				parsed+="INDD: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Do migration documenting (1) or not (0)",30) ==0) {
				parsed+="MIGD: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Do microepidemic documenting (1) or not (0)",30) ==0) {
				parsed+="MICD: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Do MLST documenting (1) or not (0).",30) ==0) {
				parsed+="MLST: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#Save sequence-allele type pairs (1) or not (0)",30) ==0) {
				parsed+="STAL: ";
				parsed+=value;
				parsed+="		";
				parsed+=temp;
				}
			else if(strncmp(temp, "#One global pop size",20) ==0) {
				cout << "---\t\tRandom population size not implemented in current version.\n"<< endl;
				}	
			else if(strncmp(temp, "#Nucleotide mismatch",10) ==0) {
				cout << "---\t\tNucleotide mismatches has been replaced with Gene-wise mismatches in current version.\n" << endl;
				}
			else if(strncmp(temp, "#Allelic mismatch",10) ==0) {
				cout << "---\t\tAllelic mismatches not implemented in current version.\n"<< endl;
				}
			else if(strncmp(temp, "#Interpopulation allelic mismatch",20) ==0) {
				cout << "---\t\tInterpopulation allelic mismatches not implemented in current version.\n"<< endl;
				}
			else if (strncmp(value, "#",0)==0) { 
				cout << "\nFirst comment line read, ending conversion." << endl;
				break;
			}
			else {
				cerr << "\nPrevious parameter line was unknown. Aborting. \n"  << endl;
				return 0;
			}

			
			

			 
		}	


	} else {
		cerr << "Error while opening input file: " << inputfn << "', " << endl; 

	}
	fclose(infile);


	cout << "\n\nFile '" << inputfn << "' converted into new format as file '" << outputfn << "'. Contents:" << endl;

	parsed+="\n\n# Output files created automatically. Adjust parameters accordingly. \n# The four letter parameter codes with the colon must be exact match. If comment ('#lorem ipsum') is too long (>180 symbols), next parameter will be read incorrectly.\n# Output probably useless if this file is set incorrectly.\n\n# Free comment space below until the three dashes (---). These and the parameters are included in the rundetails***.txt and can be used for description of the run:\n\nConverted from file '";
	parsed+=inputfn;
	parsed+="'.\n\nTo reproduce, remember to adjust seed.\n\n---\n";

	cout << "---" << endl;
	cout << parsed << endl;
	cout << "End of file contents.\n\n" << endl;




	fstream parsedfile; 
	parsedfile.open(outputfn, ios::out | ios::trunc);
	
	parsedfile << parsed << endl;

	parsedfile.close();








	return 1;

}