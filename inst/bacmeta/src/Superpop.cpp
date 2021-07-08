#include "Superpop.h"
#include "Pop.h"
#include "Para.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>


using namespace std;

vector<unordered_map<string, int>> Superpop::SeqTypes;	
unordered_map<string, int> Superpop::MLSTs;
int Superpop::allcount;
int Superpop::stcount;
int Superpop::backtoinitialmuts;
map<int,int> Superpop::indels;			
map<int,int> Superpop::recs;		
int Superpop::thisgeneration;	
	



Superpop::Superpop(int pops, int nrbacs, int nrgenes, int genelen, Para pa)  {
		
	cout << "Initializing metapopulation." << endl;

	parand=&pa; 

	cout << "\n" << endl;


	npops=pops; 
	
	ostringstream filename;
	filename   << parand->fna;
	fnadd=filename.str();
	
	nrecs=0; 
	accrecs=0; 
	ntotalmuts=0;
	ntotalaccrecs=0; 
	thisgeneration=0;
	nmigevent=0;
	
	
	if (parand->domean>0) {  
			
		fstream convfile; 
		convfile.open("outputs//MeanDistance"+ fnadd +".csv", ios::out );
		convfile << "Generation;Mean_distance" <<endl ;
		convfile << "0;0" <<endl ;
		convfile.close();
				
	}
	
	
	
	Pops.reserve(npops);
	
	stcount=0;
	allcount=0;
	

	cout << "Initializing population " <<  0 << "." << endl;

	Pop popforfor(nrbacs, nrgenes, genelen, pa); 
	Pops.push_back(popforfor);

	for (int i=0; i<npops-1; i++){	

		cout << "Initializing population " << i+1 << "." << endl;		

		Pop popforfor(i+1, nrbacs, Pops[0].bacs[0], pa); 						
		Pops.push_back(popforfor);
	}
	
	
	//Setting migration rate matrix...
	if (parand->custommigs==1) {
		migmatrix={1, vector<double>(1, 0.0)};
		migmatrix.resize(npops, vector<double>(1, 0.0));
		for (int i=0 ;i <npops; i++) {
			migmatrix[i].resize(npops,0);
		}
		
		// ...from file...
		miginit=Migmatrixupd();
	}
	else {
		for (int i=0; i<(int)Pops.size(); i++) {  // ... or from global scaling value
			Pops[i].miglambda=parand->migscaler*Pops[i].bacs.size();
		}
		
	}

	//Setting interpopulation pairwise distance computing samples.
	if (parand->popmisma>0) {  // Only global value for now. 
		distmatrix={1, vector<int>(1, 0)};
		distmatrix.resize(npops, vector<int>(1, 0));
		for (int i=0 ;i <npops; i++) {
			distmatrix[i].resize(npops,0);
		}
		
		distinit=distmatrixupd();
	}	


	SeqTypes.reserve(nrgenes);
	
	for (int i=0; i<nrgenes; i++) {
		string origen= Pops[0].parand.origenes[i];
		unordered_map<string,int> OriSeq;
		OriSeq[origen]=1;
		SeqTypes.push_back(OriSeq);
	}
	string orimlst;	
	for (int i=0; i<nrgenes; i++) {
		orimlst+="1	";
	}
	MLSTs[orimlst]=1;	
	
	
	
	for (int i=0; i<(int)Pops.size(); i++) {
			
		ostringstream filename;
		
		filename << "Pop"<<i<<"-";

		filename  << parand->fna;
		string fnadd=filename.str();
		
		if (parand->diversityinterval >0 && parand->diversityinterval <= parand->ngenerations) {
			ofstream divs;
			divs.open ("outputs//DiversityDistribution"+fnadd+".csv");
			
			divs << "Generation;Strain-id;emergeneration;snpn;Amount;Frequency" << endl;
			divs << "0;0;0;0;0;0" << endl;
			divs.close();
		}
		
		
		if (parand->mutsinterval >0 && parand->mutsinterval <= parand->ngenerations) {
			
			ofstream muts; 
			muts.open ("outputs//Mutations"+fnadd+".csv", ios::out );
			muts << "Generation;amount" <<endl ;
			muts << "0;0" <<endl ;
			muts.close();
			
		}
	
	}
	
	
	if (parand->mutsinterval >0 && parand->mutsinterval <= parand->ngenerations) {
		
		ostringstream filename;
		filename << "AllPops-";
		filename  << parand->fna;
		string fnadd=filename.str();
			
		ofstream allmuts; 
		allmuts.open ("outputs//Mutations"+fnadd+".csv", ios::out );
		allmuts << "Generation;amount" <<endl ;
		allmuts << "0;0" <<endl ;
		allmuts.close();
	}
	
	if (parand->recdocthresh > 0) {
		
		ostringstream filename;
		filename  << parand->fna;
		string fnadd=filename.str();
		
		ofstream recsfile; 
		recsfile.open ("outputs//RecEvents"+fnadd+".csv", ios::out );
		recsfile << "HD;Length;Start;Gene;Gene_ptr;Donor_ptr;Recip_ptr;Population_i;Generation" <<endl ;
		recsfile.close();
	}	
	
}

Superpop::~Superpop() { } 


void Superpop::printSuperpop(int bacinterval, int popinterval) {
	
	for (int i=0; i<npops; i++) {
		if (i % popinterval == 0) {  
			cout << "  \n   \n  -------  Population " << i << " in generation: " << thisgeneration << "  -------" << endl;
			Pops[i].printPop(bacinterval);   //Prints the genes of bacteria and some data into STDOUT. 
		}			
	}	
}


void Superpop::startmutations() {
	
	int nmuts=0; 
	for (int i=0; i<npops; i++) {
		
		nmuts=parand->getPoisrand(Pops[i].mutlambda); 
		// mutlambda=parand.mutrate*nbacs*ngenes*genelength;

		Pops[i].mutate(nmuts, thisgeneration);
		
		ntotalmuts+=nmuts; 
	}
	//cout << "Number of mutations in generation " << thisgeneration << ":   " <<  genemuts << endl;
}


void Superpop::startrecombinations() {
	
	nrecs=0;
	accrecs=0;
	for (int i=0; i<npops; i++) {

		nrecs=parand->getPoisrand(Pops[i].reclambda);
		// reclambda=parand.recrate*mutlambda;
		// mutlambda=parand.mutrate*nbacs*ngenes*genelength;

		if (parand->reclen==0) { // Do recombinations for whole genes only
			accrecs=Pops[i].recombineGene(nrecs);
		}
		else { // Do recombinations with random sized segments of genes
			accrecs=Pops[i].recombineSegment(nrecs);
		}
		
		ntotalaccrecs+=accrecs;
	}
	
}

void Superpop::startindels() {
	
	
	int nins=0; 
	int ndels=0; 
	for (int i=0; i<npops; i++) {
		
		nins=parand->getPoisrand(Pops[i].inslambda); 
		ndels=parand->getPoisrand(Pops[i].dellambda); 

		//inslambda=parand.insrate*mutlambda;
		//dellambda=parand.delrate*mutlambda;
		
		Pops[i].insertions(nins);
		Pops[i].deletions(ndels);
		
	}
	
}


bool Superpop::Migmatrixupd() { 
	// Sets new Poisson distributed migration rates to matrix where migmatrix[i][j] represents migration rate from population i to population j.

	cout << "\nCustom migration rates for population pairs enabled. Reading migration matrix from file '" << Para::migfilename << "'." << endl;
	cout << "Setting values as follows:\n" << endl;
	char filename[80];
	strcpy(filename,Para::migfilename.c_str());	
	
	ifstream infile;
	infile.open(filename);	
	vector<string> rows;
	string line;
	
	
	double temp=0;
	if (infile.is_open()) {
		for (int i =0; i<npops; i++) {
			
			getline(infile, line);
			stringstream ss(line);
			
			for (int j =0; j<npops; j++) {
				ss >> temp;	
				if (j==0 && ss.eof()) {
					cerr << "\n\nNot enough rows in '" << Para::migfilename << "' matrix. Simulation aborted." << endl; 
					return 0;
				}
				else if (ss.eof()) {
					cerr << "\n\nNot enough elements in '" << Para::migfilename << "' matrix at row " << i <<". Simulation aborted." << endl; 
					return 0;
				}
				if (i==j && temp!=0) {
					migmatrix[i][j]=0;
					printf(" %6s ", "!0");
				}
				else {
					printf(" %6.0g ", temp);
					migmatrix[i][j]=temp*Pops[i].nbacs;
				}			
			}	
			cout << endl;	
		}
	}
	else {
		cerr << "Error reading '" << Para::migfilename << "', file probably does not exist in the simulator folder. Aborting." << endl;
		return 0;
	}
	
	infile.close();
	cout << "\n" << endl;
	
	
	ifstream  inputfile(Para::migfilename, std::ios::binary);  // Document migration rate input file.
    ofstream  docfile("outputs//migrationinput"+ fnadd +".txt",   std::ios::binary);
    docfile << inputfile.rdbuf();
    
	return 1;
}	

void Superpop::migration() { 
	
	int migamount;
	for (int i=0; i<npops; i++) {
		for (int j=0; j<npops; j++) {
			if (i!=j) {
				if (parand->getZeroone()<parand->migProp) { // Check if migration happens for this route.
				
					if(parand->custommigs==0) {
						
						migamount=(parand->getPoisrand(Pops[i].miglambda)); 
						// Pops[i].miglambda=parand->migscaler*Pops[i].bacs.size();

						string temp =	Pops[i].migrate(&Pops[j], migamount);

						if (parand->migdo==1 && migamount >0) {
							ostringstream oss; 
							oss << i <<";"<< j <<";" << migamount << ";" << thisgeneration <<";"  << temp<< "\n";
							migs+= oss.str();
						}
					}	
					else {

						migamount=(parand->getPoisrand(migmatrix[i][j]));
						// migmatrix[i][j]=migrate*Pops[i].nbacs, migrate as given in migration*.input

						string temp =	Pops[i].migrate(&Pops[j], migamount);

						if (parand->migdo==1&& migamount >0) {
							ostringstream oss; 
							oss << i <<";"<< j <<";" << migamount << ";" << thisgeneration <<";"  << temp<< "\n";
							migs+= oss.str();
						}
					}	
					
					nmigevent++;
				}
			}
		}
	}
}


void Superpop::mepidemic() { 

	for (int i=0; i<(int)Pops.size(); i++) {
		int mepisperpop = parand->getPoisrand(Pops[i].bacs.size()*parand->mepiAmo);
		for (int j=0; j<mepisperpop; j++) {
			string temp = Pops[i].mepidemic();
			if (parand->micdo==1) {
				ostringstream oss; 
				oss << i <<";" << thisgeneration <<";"  << temp<< "\n";
				mics+= oss.str();
			}
		}
	} 

}


void Superpop::selection() { 

	for (int i=0; i<(int)Pops.size(); i++) {
		
		Pops[i].select();
				
	}
	
}


bool Superpop::distmatrixupd() { 

	maxdistsample=0;
	for (int i =0; i<npops-1; i++) {
		for (int j =i+1; j<npops; j++) {
			if (i==j) {
				distmatrix[i][j]=0;
				distmatrix[j][i]=0;
			}
			else {
				distmatrix[i][j]=(int)(parand->popmisma*Pops[i].nbacs+0.5);
				distmatrix[j][i]=(int)(parand->popmisma*Pops[i].nbacs+0.5);
				if(distmatrix[i][j]>maxdistsample){
					maxdistsample=distmatrix[i][j];
				}
			}			
		}		
	}

	// for (int i =0; i<npops; i++) {
	// 	for (int j =0; j<npops; j++) {
	// 		printf(" %d ", distmatrix[i][j]);
	// 	}		
	// 	cout << endl;
	// }
    
	return 1;
}	

void Superpop::newGeneration() {
	
	thisgeneration++;
	if (thisgeneration % 1000==0) { 
		cout<< "Generation: " << thisgeneration<< endl;
		
		if (parand->domean>0 && thisgeneration >= 900) {  
			
			int bacamount=0.3*(int)Pops[0].bacs.size();

			vector<int> bacsample=parand->getUtilNoRepSample(vector<int> {(int)Pops[0].bacs.size()}, bacamount);

			double meandistance=Pops[0].Burner(bacsample);
			fstream convfile; 
			convfile.open("outputs//MeanDistance"+ fnadd +".csv", ios::out | ios::app);
			convfile << thisgeneration  << ";" << meandistance  <<endl ;
			convfile.close();
			
		}
		
	}
	
	int eventorder[5] = {0,1,2,3,4};
	
	
	if (parand ->randevents) { 
		random_shuffle(eventorder, eventorder+5);
	}	
	for (int i=0; i<5; i++) {
		switch (eventorder[i]) {
			case 0:
				startmutations();
				break;
			case 1:
				startindels();
				break;
			case 2:
				startrecombinations();
				break;
			case 3:
				migration();
				break;
			case 4:
				mepidemic();
				break;
		} 
	} 
	
	document();
	
	if (parand->selerand) selection(); // If random selection of bacteria for next generation is toggled, run selection.

	
}




int Superpop::getAllele(string str, int gene) {

	allcount++;
	
	int allele;
	unordered_map<string,int>::iterator it;
	it =SeqTypes[gene].find(str);
	
	if (it != SeqTypes[gene].end()) {
		allele=it->second;
	} else {
		int newST=SeqTypes[gene].size()+1;
		SeqTypes[gene][str]=newST;
		allele=newST;
	} 
	
	return allele;	
}

int Superpop::getMLST(vector<shared_ptr<Gene>> &genes) {

	stcount++;	

	ostringstream alleles; 
	
	for (int j = 0; j< (int)genes.size(); j++) {
		alleles << genes[j]->allele <<  "	";
	}
	
	string allelestr = alleles.str();
	
	int ST;
	unordered_map<string,int>::iterator it;
	it =MLSTs.find(allelestr);
	
	if (it != MLSTs.end()) {
		ST=it->second;
	} else {
		ST=MLSTs.size()+1;
		MLSTs[allelestr]=ST;
	}

	return ST;	
}



void Superpop::document() { 
	

	if (parand->genewiseval >0 && thisgeneration % parand->genewiseval == 0  && parand->popeval >0 && thisgeneration % parand->popeval == 0) { 
		

		vector<vector<int>> popsamples;
		popsamples.resize(npops, vector<int>(1, -1));

		for (int i=0 ;i <npops; i++) {
			popsamples[i].resize(maxdistsample,-1);
		}

		for (int i=0; i<(int)Pops.size(); i++) {

			int intersize=popsamples[i].size();
			int intrasize=parand->genewisemisma*(int)Pops[i].bacs.size();

			if (intrasize>intersize) {

				vector<int> intrasample=parand->getUtilNoRepSample(vector<int> {(int)Pops[i].bacs.size()}, intrasize);
				Pops[i].GeneWiseDiffs(intrasample);

				vector<int> intersample=parand->getUtilNoRepSample(intrasample, intersize);
				popsamples[i]=intersample;

			}
			else if (intrasize<intersize) {
				vector<int> intersample=parand->getUtilNoRepSample(vector<int> {(int)Pops[i].bacs.size()}, intersize);
				popsamples[i]=intersample;

				vector<int> intrasample=parand->getUtilNoRepSample(intersample, intrasize);
				Pops[i].GeneWiseDiffs(intrasample);

			}
			else {
				vector<int> samesample=parand->getUtilNoRepSample(vector<int> {(int)Pops[i].bacs.size()}, intersize);
				popsamples[i]=samesample;

				Pops[i].GeneWiseDiffs(samesample);
			}

		}



		for (int i=0; i<(int)Pops.size()-1; i++) {
			for (int j=i+1; j<(int)Pops.size(); j++) {
				int samplesize=distmatrix[i][j];
				if (samplesize>0) {
					Pops[i].InterpopGeneWiseDiffs(popsamples[i], popsamples[j], Pops[j]);
				}
				// printf("%d,%d: %d ", i, j, samplesize);
			}		
		// cout << endl;
		}


	}
	else if (parand->genewiseval >0 && thisgeneration % parand->genewiseval == 0) { 
		for (int i=0; i<(int)Pops.size(); i++) {

				int bacamount=parand->genewisemisma*(int)Pops[i].bacs.size();
				vector<int> bacsample=parand->getUtilNoRepSample(vector<int> {(int)Pops[i].bacs.size()}, bacamount);

				Pops[i].GeneWiseDiffs(bacsample);
		}
	}
	else if (parand->popeval >0 && thisgeneration % parand->popeval == 0) { 
		
		vector<vector<int>> popsamples;
		popsamples.resize(npops, vector<int>(1, -1));
		for (int i=0 ;i <npops; i++) {
			popsamples[i].resize(maxdistsample,-1);
		}

		for (int i=0; i<(int)Pops.size(); i++) {

			int bacamount=popsamples[i].size();
			vector<int> bacsample=parand->getUtilNoRepSample(vector<int> {(int)Pops[i].bacs.size()}, bacamount);
			popsamples[i]=bacsample;

		}



		for (int i=0; i<(int)Pops.size()-1; i++) {
			for (int j=i+1; j<(int)Pops.size(); j++) {
				int samplesize=distmatrix[i][j];
				if (samplesize>0) {
					Pops[i].InterpopGeneWiseDiffs(popsamples[i], popsamples[j], Pops[j]);
				}
				// printf("%d,%d: %d ", i, j, samplesize);
			}		
		// cout << endl;
		}

		
	}


	
	if  (parand->debugval>0 && thisgeneration % parand->debugval==0) { 
		int popi=parand->debugpop;
		cout << parand->debugval<<  "generation , pop" <<  popi << endl;	
		Pops[popi].PopInfo();
	}
				
	
	if (parand->summaryinterval >0 && thisgeneration % parand->summaryinterval == 0) { 
		for (int i=0; i<(int)Pops.size(); i++) {
			Pops[i].summary(thisgeneration, i);
		}	
	}
	
	if (parand->diversityinterval >0 && thisgeneration % parand->diversityinterval == 0) { 
		
		for (int i=0; i<(int)Pops.size(); i++) {
			Pops[i].Diversity(thisgeneration, i);
		}	
	}
	
	if (parand->recsinterval >0 && thisgeneration % parand->recsinterval == 0) { 
		
			
			ostringstream filename;
			
			filename <<thisgeneration<<"-";

			filename  << parand->fna;
			string fnadd=filename.str();
			
			
			ofstream recsfile; 
			recsfile.open ("outputs//RecCounts"+fnadd+".csv");
			
			recsfile << "Effect;Amount" << "\n"; 	

			for (auto const &recc : recs) {
				recsfile << recc.first<< ";" << recc.second << "\n"; 	
			}
			
			recsfile.close();  
		
	}
	
	if (parand->mutsinterval >0 && thisgeneration % parand->mutsinterval == 0) { 
		
		int allpopmuts=0;
		
		for (int i=0; i<(int)Pops.size(); i++) {
			
				
			ostringstream filename;
			
			filename << "Pop"<<i<<"-";

			filename  << parand->fna;
			string fnadd=filename.str();
			
			ofstream mu; 
			mu.open ("outputs//Mutations"+fnadd+".csv", ios::out | ios::app);
			
			mu << thisgeneration  << ";" <<  Pops[i].tmuts <<endl ;

			mu.close();
			
			allpopmuts+=Pops[i].tmuts;
			
		}	
		
		ostringstream filename;
			
		filename << "AllPops"<<"-";

		filename  << parand->fna;
		string fnadd=filename.str();
		
		ofstream mu; 
		mu.open ("outputs//Mutations"+fnadd+".csv", ios::out | ios::app);
		
		mu << thisgeneration  << ";" <<  allpopmuts <<endl ;

		mu.close();  
		
	}
	
	if (parand->seqval >0 && thisgeneration % parand->seqval == 0) {  
		ofstream seqtofile; 
		ostringstream filename;
		
		filename << "outputs//" << "Sequences" << thisgeneration << "-" << parand->fna  <<".csv";
		
		string tofile=filename.str();
		
		seqtofile.open (tofile);
		seqtofile <<  "Population;Amount";
		
		
		double proportion=parand->seqsample;
		
		for (int p=0; p<(int)Pops.size(); p++) {
				
			int bacamount=proportion*(int)Pops[p].bacs.size();
		
			vector<int> randbacinds=parand->getUtilNoRepSample(vector<int> {(int)Pops[p].bacs.size()}, bacamount);
			
				
			
			map<shared_ptr<Bac>,int> bseqs;
			
			for (auto & v : randbacinds) {
				
				auto it = begin(bseqs);
				 
				it = bseqs.find(Pops[p].bacs[v]);
				
				if (it!=bseqs.end()) {
					
					it->second=it->second+1;
					
				}
				else {
					bseqs[Pops[p].bacs[v]]=1;
				}
				
			
			
			}	
		
			for (int i=0; i<parand->ngenes; i++) {
				seqtofile <<  ";";
				seqtofile <<  "Gene" << i;
			}
			seqtofile << endl;
			
			for (auto & t: bseqs) {
				seqtofile << p << ";" << t.second << ";";
				for (int l=0; l<parand->ngenes; l++) {
					seqtofile << t.first->genes[l]->getGeneText(1) << ";" ;	
				}
				seqtofile <<  endl;
			}
			
		}
		seqtofile.close(); 
	}
}

void Superpop::savedocumentation() { 
	
	for (int i=0; i<(int)Pops.size(); i++) {
		
		
		if (parand->summaryinterval >0 && parand->summaryinterval <= parand->ngenerations) {
			ofstream summarytofile; 
			summarytofile.open ("outputs//Summaries"+fnadd+".txt");
			if (parand->mlstdo) {
				summarytofile << "Number of observed genotypes: " << MLSTs.size() << " \n" << endl;	
			}
			for (int k=0; k<(int)Pops.size(); k++) {
					summarytofile << Pops[k].summarystr << endl;	
					summarytofile << "\n  \n   =================================== \n " <<endl;
			}
			
			summarytofile << "Total number of migration events " << nmigevent << "  \n"	<< endl;
			summarytofile.close(); 
		}
		
	}

	if (parand->indeldo==1) {
		ofstream inds; 
		inds.open ("outputs//Indels"+fnadd+".csv");
		inds << "length;amount"<< endl;
		for (auto const &ent : indels) {
			inds << ent.first<< "; " << ent.second << "\n"; 	
		}
		
		inds.close(); 
	}
	
	if (parand->micdo==1) {
		ofstream micros; 
		micros.open ("outputs//Microepidemics"+fnadd+".csv");

		micros << "Population;Generation;MLST;Size" << " \n" ;
		micros << mics << " \n" ;
		micros.close(); 
	}
	
	if (parand->migdo==1) {
		ofstream migras; 
		migras.open ("outputs//Migrations"+fnadd+".csv");
		migras <<  "SourcePopulation;" <<"DestinationPopulation;" << "MigrationAmount;" << "Generation;" << "ST_1-amount,ST_2-amount..." << endl;
		migras << migs << " \n" ;
		migras.close(); 
	}
	if (parand->mlstdo>0) {
		ofstream MLSTtofile; 
		MLSTtofile.open ("outputs//MLST"+fnadd+".txt");
		for (int k=0; k<(int)Pops.size(); k++) {
			for (int i=0; i < (int)Pops[k].bacs.size() ; i++) {
				MLSTtofile << Pops[k].bacs[i]->MLST << "	" ;
				for (int j = 0; j< (int)Pops[k].bacs[i]->genes.size(); j++) {
					MLSTtofile << Pops[k].bacs[i]->genes[j]->allele <<  "	";
				}
				MLSTtofile <<" \n" ;
			}
		}
		MLSTtofile.close(); 
		
		
		if (parand->alleseq>0) {
			
			
			 vector<map<int, string>> stvec;
			stvec.reserve((int)SeqTypes.size());
			for (int i=0; i<(int)SeqTypes.size(); i++) {
				map<int, string> temp;
				for (unordered_map<string,int>::iterator it=SeqTypes[i].begin(); it!=SeqTypes[i].end(); ++it) {
					temp[it->second]= it->first; 
					//STtofile << i << " " << it->second <<  " " << it->first ;
					//STtofile <<" \n" ;		
				}
				stvec.push_back(temp);
			}
			
			
			ofstream ordSTtofile; 
			ordSTtofile.open ("outputs//MLSTsequences"+fnadd+".txt");
			for (int i=0; i<(int)stvec.size(); i++) {
				
					ordSTtofile << "Gene: " << i << endl;
				for (auto& x: stvec[i]) {
					ordSTtofile << x.first <<  "	" << x.second << " \n";
				}
				
					ordSTtofile << "\n" << endl;
			}
			ordSTtofile.close();  
		}
	
	}
	
}
