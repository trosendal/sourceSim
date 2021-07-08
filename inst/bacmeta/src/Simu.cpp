#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <random>
#include "Gene.h"
#include "Bac.h"
#include "Para.h"
#include "Pop.h"
#include "Superpop.h"
#include <time.h> 
#include <algorithm>
#include <map>

using namespace std;
/* Bacmeta -Simulator of genetic structure of bacteria strain in multiple populations 

To run:

1. Compile with attached makefile, non-standard libraries not required
2. Set parameters for the simulation in "simu.input" -file
3. Run with "./simu"  

Short documentation available in "README.md"
*/

string Para::inputfilename;
string Para::outmod;
string Para::migfilename;


int argparser(int argc, char* argv[], string& inputfile, string& outputfile, string& migfile, string& flagsused) {  //Handle optional command line flag arguments. 
	
    if (argc < 2) {
        cout << "\n Using default simu.input-file for parameter inputs. \n\n" << endl;
    }
    else if (argc < 3) {
		cout << "\nArguments are to be given with flag, e.g. \"./simu -p 123\" to determine specific parameter input file simu123.input to be used." << endl;
		if (string(argv[1]) == "-h") {
				cout << "\nHelp for arguments: \n" << endl;
				cout << "-p   Specify 'simu---.input' file for parameter input.  The argument replaces '---'. \nFor example './simu -i 123' would use simu123.input.\n"<< endl;
				cout << "-m   Specify 'migration---.input' file for  migration rate matrix input.  The argument replaces '---'. \nFor example './simu -m 123' would use migration123.input.\n" << endl;
				cout << "-b   Specify both of the files above for parameter input.  The argument replaces '---'. \nFor example './simu -b 123' would use simu123.input AND migration123.input. \n\n" << endl;
				cout << "-o   Specify modifier for the output filenames. \nFor example './simu -p 123 -o 456' would use simu123.input BUT output files would have filename modifier 456, ignoring the value in simu123.input. \n\n" << endl;
				cout << "Simulator will default to 'simu.input' and/or to 'migration.input', if no files are specified with these flags." << endl;
				
				cout << "Arguments work at least with numbers and alphabetic letters. Symbols have NOT been tested thoroughly enough to be advised." << endl;
				
				return 0;
			}
	}
	else {
		for (int i = 1; i < argc; ++i) {
			
			if (string(argv[i]) == "-p") {
				if (i + 1 < argc) { 
					i++;
					inputfile = "simu"+string(argv[i])+".input"; 
					
					struct stat inputinfo;

					if( stat( inputfile.c_str(), &inputinfo ) != 0 ) {
						cerr << "No input file '" << inputfile << "' found. Aborting. \n" <<endl;
						return 0;
					}
					else  {
						cout << "Using specified '" << inputfile << "' for parameter input. \n\n" << endl;
					}
	
				} else { 
					  cerr << "Argument missing for parameter input file. Aborting" << endl;
					  return 0;
				}  
			}
			 
			else if (string(argv[i]) == "-m"){
				if (i + 1 < argc) { 
					i++;
					migfile = "migration"+string(argv[i])+".input"; 
					
					struct stat inputinfo;

					if( stat( migfile.c_str(), &inputinfo ) != 0 ) {
						cerr << "No migration input file '" << migfile << "' found. Aborting. \n" <<endl;
						return 0;
					}
					else  {
						cout << "Using specified '" << migfile << "' for migration matrix input. \n\n" << endl;
					}
	
				} else { 
					  cerr << "Argument missing for migration file. Aborting" << endl;
					  return 0;
				}  
			}
			
			else if (string(argv[i]) == "-b"){
				if (i + 1 < argc) { 
					i++;
					inputfile = "simu"+string(argv[i])+".input"; 
					migfile = "migration"+string(argv[i])+".input"; 
					
					
					struct stat inputinfo;

					if( stat( inputfile.c_str(), &inputinfo ) != 0 ) {
						cerr << "No input file '" << inputfile << "' found. Aborting. \n" <<endl;
						return 0;
					}
					else  {
						cout << "Using specified '" << inputfile << "' for parameter input. \n\n" << endl;
					}
					
					struct stat inputinfo2;
					
					if( stat( migfile.c_str(), &inputinfo2 ) != 0 ) {
						cerr << "No migration input file '" << migfile << "' found. Aborting. \n" <<endl;
						return 0;
					}
					else  {
						cout << "Using specified '" << migfile << "' for migration matrix input. \n\n" << endl;
					}
	
				} else { 
					  cerr << "Argument missing for parameter input file. Aborting" << endl;
					  return 0;
				}  
			}
			
			else if (string(argv[i]) == "-o"){
				if (i + 1 < argc) { 
					i++;
					
					outputfile =  string(argv[i]);
					
					cout << "Using output file name modifier: '" << outputfile << "'. \n\n" << endl;
	
				} else { 
					  cerr << "Argument missing for output file. Aborting" << endl;
					  return 0;
				}  
			}
			else if (string(argv[i]) == "--convert"){  // This flag exists only for legacy reasons, and should be used only if you want to convert old style input files into new format and know the process. Contact author for help if such conversion is necessary!
				if (i + 2 < argc) { 

					i++;
					char inputfn[80]; 
					strcpy (inputfn,  string(argv[i]).c_str());

					i++;
					string outputfn = string(argv[i]);

					cout << "\nConverting old style input file: '" << inputfn << "' into new input format as while: '"<< outputfn << "'. \n\n" << endl;

					if (Para::parser(inputfn, outputfn)!=1) {

						cout << "Input file conversion failed." << endl;

						return 2;
					}
					
					cout << "Input file conversion done." << endl;

					return 2;

	
				} else { 
					  cerr << "Argument missing for converting, give '--convert inputfilename outputfilename'. Aborting" << endl;
					  return 0;
				}  
			}
			
			
			
			else {
				cerr << "Unknown " << argv[i] << " option. Aborting" << endl;
				return 0;
			}
		}

		for (int i = 1; i < argc; ++i) {
			flagsused+=string(argv[i])+" ";
		}
	}	
	
	return 1;
}





//
// MAIN 
//


int main(int argc, char* argv[])
{
	clock_t tim;  // Simple run time documentation. 
	tim= clock();
	
	
	string inputfile="simu.input";  //Variables for parameter input filenames...
	string outputfile="";
	string migfile="migration.input";
	string flagsused="";


	int argparserresult=argparser(argc, argv, inputfile, outputfile, migfile, flagsused);
	if (argparserresult==0) { //... updated if command line arguments are used. 
		return 0;  // Aborting if argparser returns error because of input file problem.
	} else if (argparserresult==2){ 
		return 0;
	}
	
	//Parameter input files are always of the form simu*.input. 
	//Migration input files are always of the form migration*.input
	//Outputfile modifier is alphanumeric string attached to all outputfiles to differentiate between runs.
	
	
	
	
	// Creation of Para object for parameter input and random number generation. More info in the class itself.
	Para::inputfilename=inputfile;
	Para::outmod=outputfile;
	Para::migfilename=migfile; 
	
	Para mainparand; 

	if (Para::outmod=="error") {
		return 0;
	}
	

	mainparand.initZipf(1.0);	 
	// The parameters can now be requested from the Para object in scope.
	
	

	struct stat info;  // Create ouput folder if needed.
	if( stat( ".//outputs", &info ) != 0 ) {
		printf( "No 'outputs' folder found. Creating one.\n" );
		mkdir("outputs", S_IRWXU);
		printf( "'outputs' folder created\n" );
	}
	
	
	cout<< "\nReading and initializing successful:" << endl;  
	cout<< "Simulation run file name modifier: "<< mainparand.fna << endl;
	cout << "Number of generations: " << mainparand.ngenerations << endl;
	cout << "Seed parameter: " <<  Para::trueseed << " \n" << endl;
	
	
	
	string runtofile="";  //Save simu*.input file contents, so next run can modify it, even before documentation for previous run takes place.
	string line;
	ifstream inputs (inputfile);
	if (inputs.is_open()) {
		while ( getline (inputs,line)){
			if (line.substr(0,1)=="-")	{
				runtofile  += line;
				runtofile  += '\n';
				break;
			}
			
			else {
				runtofile  += line;
				runtofile  += '\n';
			}
		}
	}
    inputs.close();
	
	
	
	// Creation of Superpop i.e. the metapopulation object that handles the collection of all the populations in the simulation. 
	Superpop superpop(mainparand.npops, mainparand.nbacs, mainparand.ngenes, mainparand.genelen, mainparand);


	if (mainparand.custommigs && !superpop.miginit) { // Check that initializing of migration rates from migration*.input was succesful if requested, and abort run if it failed.
		return 0;
	}

	//Simulation itself happens here
	for (int i=0; i<mainparand.ngenerations ; i++) {
		superpop.newGeneration(); 
	}
	
	
	superpop.savedocumentation();// To run the requested documentations of simulation.
	
	
	
	//Document the run so it can be repeated.
	ostringstream filename;  // Trickery for portable output filename handling. 
	filename  << mainparand.fna;
	string fnadd=filename.str();
	
	fstream runfile; 
	runfile.open("outputs//rundetails"+fnadd+".txt", ios::out | ios::trunc);
	runfile << runtofile; // runtofile contains the contents of simu*.input used.
	
	//Include the used seed, especially useful when using random seed from time. Also simple runtime documentation.
	runfile << "\n\n#Details of simulation run: \n" << endl;
	runfile<< "\n#Seed used: " << Para::trueseed;
	if (flagsused!=""){
		runfile << " \n#Flags used: " << flagsused << endl;
	}
	else {
		runfile << " \n#Flags used: none " << endl;
	}
	runfile << " \n#Run time of simulation in approximate real time:   " << float((clock() -tim))/CLOCKS_PER_SEC << "s" << endl;


	runfile.close();
	
	
	// Done
	cout << "\nSimulation "<< mainparand.fna <<" finished and output files created." << endl;
	

 	return 0;
}


