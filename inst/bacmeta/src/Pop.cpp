#include "Pop.h"
#include "Superpop.h"
#include "Bac.h"
#include "Para.h"
#include <random>
#include <algorithm>
#include <map>
#include <forward_list>
#include <sstream>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <stdint.h>

using namespace std;


Pop::Pop(int nbacteria, int nrgenes, int genelen, Para pa)  {
	popindex=0;

	shared_ptr<Bac> FirstBac = make_shared<Bac>(nrgenes, genelen, pa);  // Creating new bacterium with randomly generated sequences as loci.
	//shared_ptr's are used for the superior efficiency and usability.

	parand=pa; // Setting Para object parand for the use of whole class


	nbacs=nbacteria;  
	ngenes=nrgenes;
	genelength=genelen;
	mutlambda=parand.mutrate*nbacs*ngenes*genelength;
	reclambda=parand.recrate*mutlambda;
	inslambda=parand.insrate*mutlambda;
	dellambda=parand.delrate*mutlambda;
	
	nrecombinations=0; 
	nrectries=0; 
	tmuts=0;
	snpsfromuts=0;
	mutsonsnp=0; 
	mutstoanc=0; 
	nullrecombs=0;
	recsnps=0;
	dsnps=0;
	addrecn=0;
	delrecn=0;
	bacs.reserve(nbacs);
	GCs=FirstBac->GC;  
	aa=0;  
	cc=0;		
	gg=0;	
	tt=0;
	bacsimmigrated=0;  
	bacsemigrated=0;  
	nmepidems=0;
	nmpedimbacs=0;
	mismastr="";
	allemismastr="";


	for (int i=0; i<nbacs; i++){
		bacs.push_back(FirstBac);  // Add rest of bacteria as copies from first bacterium bacsforfor  
	}
	
	for (int j=0; j<ngenes; j++) {
		Para::origenes.push_back(bacs[0]->genes[j]->genstr);
	}	
	
	ostringstream filename;
	filename  << pa.fna;
	string fnadd=filename.str();
	
	
	if (parand.iniseq==1) {
		ofstream ori; 
		ori.open ("outputs//InitialSequences"+fnadd+".txt");
		
		for (int i=0; i < nrgenes ; i++) {
			ori << parand.origenes[i] << " \n  \n"  ; // write initial genes to file
		}
		ori.close();  
	}
	
}

// Constructs population based on the bacterium object (Bac) given as argument.
Pop::Pop(int popi, int nbacteria, shared_ptr<Bac> bacinput, Para pa)  { 

	popindex=popi;

	parand=pa;
	nbacs=nbacteria;
	ngenes=bacinput->ngenes;
	genelength=bacinput->genelength;
	mutlambda=parand.mutrate*nbacs*ngenes*genelength;
	reclambda=parand.recrate*mutlambda;
	inslambda=parand.insrate*mutlambda;
	dellambda=parand.delrate*mutlambda;
	
	GCs=bacinput->GC;
	tmuts=0;
	snpsfromuts=0;
	mutsonsnp=0; 
	mutstoanc=0; 
	bacsimmigrated=0;  
	bacsemigrated=0;  
	nmepidems=0;
	nrecombinations=0;
	nrectries=0; 
	nmpedimbacs=0;
	nullrecombs=0;
	recsnps=0;
	dsnps=0;
	addrecn=0;
	delrecn=0;
	mismastr="";
	allemismastr="";
	 
	bacs.reserve(nbacs);
	
	shared_ptr<Bac> FirstBac = make_shared<Bac>(*bacinput);
	for (int i=0; i<nbacs; i++) {
		bacs.push_back(FirstBac);
	}
	
}


void Pop::printPop(int interval) {  
	
	cout << "Amount of bacteria migrated in: " << bacsimmigrated << endl;	
	for (int i=0; i<(int)bacs.size(); i++) {
		if (i % interval==0) { 
			cout << "  \n        Genes of bacterium " << i << endl;
			bacs[i]->printGenes();
		}
		
	}		
}

void Pop::mutate(int nmuts, int thisgeneration) {  
	int mutbac=0;
	tmuts+=nmuts;
	int snpd;

	for (int i=0; i<nmuts; i++) {	
		mutbac = bacs.size()*parand.getZeroone();
		
		//~ cout << "Before mutation \n\n";
		//~ BacRangeInfo(mutbac, mutbac+1);
		
		shared_ptr<Bac> MutBac = make_shared<Bac>(*bacs[mutbac]);
		MutBac->mutate(parand, snpd);
		MutBac->emergeneration=thisgeneration;
		bacs[mutbac]=MutBac;
		
		//~ cout << "After mutation \n\n";
		//~ BacRangeInfo(mutbac, mutbac+1);
		if (snpd==1) {
			snpsfromuts++;
		}
		else if (snpd==0) {
			mutsonsnp++;
		}
		else {
			mutstoanc++;
		}			
	}	
	
	
	
		
}

int Pop::recombineSegment(int nrecs) {  // Segment of loci recombination
	
	//Variables for:
	int recbac=0; // index of recipient bacterium
	int donbac=0; //index of donor bacterium
	int recgene=0; // index of recombination gene
	int rocus=0; // Starting point of recombination
	int reclength=0; // Length of recombination sequence -1 !! so that rocus + reclength is the end of recombinations, not one over
	int accrecs=0; // Amount of successful recombination tries.
	nrectries+=nrecs;
	
	vector<pair< int, string>>::iterator it;
	map<int, string>::iterator donlowi,donhighi,reclowi,rechighi;
	vector<pair< int, string>> difs;
	int inshd;
	
	vector<int>::iterator dit;
	set<int>::iterator donlowd,donhighd,reclowd,rechighd;
	vector<int> deldifs;
	
	int delhd;
	
	
	map<int, char>::iterator donlows,donhighs,reclows,rechighs;
	map< int, char>::iterator sit;
	map< int, char>::iterator sit2;
	
	int snphd;
									
	
	for (int i=0; i<nrecs; i++) {		
		donbac = bacs.size()*parand.getZeroone();
		recbac = bacs.size()*parand.getZeroone();
		recgene = ngenes*parand.getZeroone();
		
		if (bacs[donbac]->genes[recgene].get()==bacs[recbac]->genes[recgene].get()) { // Skip recombination actions if same allele
			nullrecombs++;
			Superpop::recs[0]++;
		} 
		else { 
		
			int totalsnps=bacs[donbac]->genes[recgene]->snps.size()+bacs[recbac]->genes[recgene]->snps.size();
			int totalindels=bacs[donbac]->genes[recgene]->dels.size()+bacs[recbac]->genes[recgene]->dels.size() + bacs[donbac]->genes[recgene]->inslens+bacs[recbac]->genes[recgene]->inslens;
			
			if (totalsnps==0 && totalindels==0) { // Skip recombination actions if original alleles
				nullrecombs++;
				Superpop::recs[0]++;
			}
			
			else {		// "Proper" recombination now with indel-extravaganza, hold your hat.
				
				double tryrec=parand.getZeroone();
				double hd=0; 
				
				int lenmod=0;
				
				reclength = parand.getGeorand();
				
				rocus = (genelength+reclength)*parand.getZeroone()-reclength; // Trickery to include recombinations starting before the gene
				
				if( rocus < 0 ) {
					lenmod=rocus;
					rocus = 0;
				}
				
					
				char oldbase = bacs[donbac]->genes[recgene]->genstr[rocus];  // Get recombination rocus i.e. starting point thats not on deletion-site
				
				while (oldbase=='-') { 
					if (rocus ==0) { // Beginning of gene deleted
						break;
					}
					else { // Starting point on deletion-site -> get new
						rocus = (genelength+reclength)*parand.getZeroone()-reclength;
						if( rocus < 0 ) {
								lenmod=rocus;
								rocus = 0;
						}
					}
					oldbase = bacs[donbac]->genes[recgene]->genstr[rocus];
				}
					
				reclength += lenmod;  // Get the length of overlap in the gene
			

				string donseq=bacs[donbac]->genes[recgene]->genstr.substr(rocus, reclength); 
								
				while (donseq[donseq.length()-1]=='-') { // ...remove excess deletion sites from recsite end...
					donseq.erase(donseq.end()-1);
				}	

				
				reclength=donseq.length() ; // and get the actual reclength

				
				int recend=rocus+reclength -1; // Take recombination endpoint
				
				
				string recipientseq = bacs[recbac]->genes[recgene]->genstr.substr(rocus, reclength); // String of sequence to be replaced
				
				// BEGIN Hamming Distance COMPUTING
				
				
				// hd from SNPS
				
				donlows=bacs[donbac]->genes[recgene]->snps.lower_bound(rocus);  // Insertions in the donor sequence
				donhighs=bacs[donbac]->genes[recgene]->snps.upper_bound(recend);
				
				reclows=bacs[recbac]->genes[recgene]->snps.lower_bound(rocus); // Insertions in the recipient sequence
				rechighs=bacs[recbac]->genes[recgene]->snps.upper_bound(recend);
				
				sit=donlows;
				sit2=reclows;
				
				
				snphd=distance(donlows,donhighs)+distance(reclows,rechighs);
				
				
				while (sit != donhighs && sit2 !=rechighs){
					if (sit->first < sit2->first) {
						++sit;
					} else if (sit2->first < sit->first) {
						++sit2;
					} else { // equal keys
						if (sit->second==sit2->second) {
							snphd-=2;
						}
						if (sit != donhighs) {
							++sit;
						}
						if (sit2 !=rechighs){
							++sit2;
						}
					}
					
					if (sit == donhighs) {
						while (sit2 !=rechighs){
							if (sit->first == sit2->first && sit->second==sit2->second) {
							snphd-=2;
							}
							++sit2;
						}
					} else if (sit2 ==rechighs){
						while (sit !=donhighs){
							if (sit->first == sit2->first && sit->second==sit2->second) {
							snphd-=2;
							}
							++sit;
						}
					} 
					
				}
				
				hd+=snphd;  
				

				inshd=0;
				delhd=0;
				
				if (totalindels >0 ){
					
					if (tryrec > (double)(pow(10,-(hd/(donseq.length()))*parand.recacc))) {
						Superpop::recs[-1]++;	
						goto nextrec;
					}
					// first hd from insertions
					
					// Get iterators for the insertion ranges that are relevant here
					donlowi=bacs[donbac]->genes[recgene]->ins.lower_bound(rocus);  // Insertions in the donor sequence
					donhighi=bacs[donbac]->genes[recgene]->ins.upper_bound(recend);
					
					reclowi=bacs[recbac]->genes[recgene]->ins.lower_bound(rocus); // Insertions in the recipient sequence
					rechighi=bacs[recbac]->genes[recgene]->ins.upper_bound(recend);
					
					difs.resize(distance(donlowi, donhighi)+distance(reclowi,rechighi)); // Vector for differing insertions
					
					it=set_symmetric_difference(donlowi,donhighi,reclowi,rechighi, difs.begin()); // Get differing insertions, it points to last differing ins, i.e. last non null difs element.
					
					
					if (it!=difs.begin()) {  // If there is differences...
						difs.resize(it-difs.begin());			// drop null elements
					
						bool lastunique=1;
					
						for (int j=0; j<(int)difs.size()-1; j++ ) {
							
							if (difs[j].first==difs[j+1].first) { // if there is two inertions with same starting point (ins[n].first) compare.
								
								// Go through the string until either ends...
								for (int k= 0; k<(int)min(difs[j].second.length(),difs[j+1].second.length()); k++) {
									if (difs[j].second[k] != difs[j+1].second[k]) {
										inshd++;
									}
									
								}
								inshd+=fabs(difs[j].second.length()-difs[j+1].second.length()); // ...and add the diffence in length
								
								if (j+1==(int)difs.size()-1){
										lastunique=0;
								}
								
								j++;
							}
							
							else { // No insertion at same start point, add length of insertion to hd.
								inshd+=difs[j].second.length();
							}

							// If deemed preferable in the future, implementation could be added to reduce length of recombination in original sequence as more insertions are donated. 
							//~ if (difs[j].first <rocus+reclength-inshd) { 
								//~ break;
							//~ }
						}
						if (lastunique) {
							inshd+=difs[difs.size()-1].second.length();
						}
					}
					
					hd+=inshd;
					
					if (tryrec > (double)(pow(10,-(hd/(donseq.length()))*parand.recacc))) {
						Superpop::recs[-1]++;	
						goto nextrec;
					}
				
					// hd from deletion
					
					// Get iterators for the insertion ranges that are relevant here
					
					donlowd=bacs[donbac]->genes[recgene]->dels.lower_bound(rocus);  // deletions in the donor sequence
					donhighd=bacs[donbac]->genes[recgene]->dels.upper_bound(recend);
					
					reclowd=bacs[recbac]->genes[recgene]->dels.lower_bound(rocus); // deletions in the recipient sequence
					rechighd=bacs[recbac]->genes[recgene]->dels.upper_bound(recend);
					
					
					deldifs.resize(bacs[donbac]->genes[recgene]->dels.size()+bacs[recbac]->genes[recgene]->dels.size()); // Vector for differing deletions
					
					dit=set_symmetric_difference(donlowd, donhighd, reclowd,rechighd, deldifs.begin()); // Get differing deletions, dit points to last differing dels 
					delhd=distance(deldifs.begin(),dit);
					
					hd+=delhd;
					
					
				}	
				
				// END OF Hamming Distance COMPUTING
				
				hd=(double)(hd/(donseq.length()));  // Adjust with length
							
				if (hd==0) { // Skip recombination actions if original alleles
					nullrecombs++;
					Superpop::recs[0]++;
				}
				
				else if (tryrec<(double)(pow(10,-hd*parand.recacc))) { 
					//Partial gene recombination
													  
												  // if (delhd> 0) {
													  
													 //  cout << "\n\n\n Recombination  \n";
													
														
													 //  cout << "hd: " << hd*donseq.length() << "   rel hd:" <<  hd << "  length:" << donseq.length() << endl;
													 //  cout << "gene: " << recgene << "  site: " <<  rocus << "-" << recend << endl;

													 //  cout << "\n\n\n     donor reci";
													 //  BacRangeInfo(donbac, donbac, recgene);
													
													 //  cout << "\nBefore rec reci\n\n";
													 //  BacRangeInfo(recbac, recbac, recgene);
														
												  // }
						  
						
						
						if (parand.recdocthresh>0 && parand.recdocthresh <= hd*(donseq.length())) {
						
							ostringstream filename;
							filename  << parand.fna;
							string fnadd=filename.str();
							
							ofstream recsfile; 
							recsfile.open ("outputs//RecEvents"+fnadd+".csv", ios::out | ios::app);
							recsfile << hd*(donseq.length()) <<";"<< (donseq.length()) << ";" << rocus << ";" << recgene << ";" <<bacs[donbac]->genes[recgene] << ";" << bacs[donbac] << ";" << bacs[recbac]<<";"<< popindex <<";"<< Superpop::thisgeneration <<endl ;
							recsfile.close();
							
							
						}
						
						
						
						shared_ptr<Gene> RecGene = make_shared<Gene>(*bacs[recbac]->genes[recgene]);
												
						int snpcountchange=0-(RecGene->snps.size());
						
						// TRANSFER INSERTS 
						
						if (inshd!=0) {
							
							reclowi=RecGene->ins.lower_bound(rocus); // Get iterators to the new shared_ptr thats created for the new allele
							rechighi=RecGene->ins.upper_bound(recend);
						
							if (0<distance(reclowi,rechighi)) {  		// And erase insertions in between
								RecGene->ins.erase(reclowi,rechighi);  // inslens adjustment done in the Bac::recombine
							}
							if (0<distance(donlowi,donhighi)) {      // Add donated ones. 
								RecGene->ins.insert(donlowi,donhighi);
							}
							
							if (parand.recsites==1) { // Delete previous recsites at insert sites. 
												// Bac::recombine handles inserting new ones 
								if (RecGene->insrecs.size()>0) {  
									
									auto recit = begin(RecGene->insrecs);
				
									while (recit != end(RecGene->insrecs)) {
										
										if (rocus+reclength < recit->first ) {
											break;
										}
										else if (rocus < recit->first )
											recit = RecGene->insrecs.erase(recit);
										else
											++recit;
									}
								}
							}
						}
						
						// TRANSFER DELETIONS		
						
							
					
						if (delhd!=0) {
							
							reclowd=RecGene->dels.lower_bound(rocus); // Get iterators to the new shared_ptr thats created for the new allele
							rechighd=RecGene->dels.upper_bound(recend);
						
							if (0<distance(reclowd,rechighd)) {  		// And erase insertions in between
								RecGene->dels.erase(reclowd,rechighd);  // inslens adjustment done in the Bac::recombine
							}
							if (0<distance(donlowd,donhighd)) {      // Add donated ones. 
								RecGene->dels.insert(donlowd,donhighd);
							}
							
						}
						
									
						
						// TRANSFER SNPS
						if (snphd!=0) { 
							map<int,char>::iterator itlow,itup;			
							
							// Delete snps of recipient from the recombination site
							itlow=RecGene->snps.lower_bound(rocus);  // Get range of snp deletions.
							itup=RecGene->snps.upper_bound(recend); 
								
							RecGene->snps.erase(itlow,itup);        // erases all the snps that are in the range.
								
							RecGene->snps.insert(donlows,donhighs);
							
							snpcountchange+=RecGene->snps.size();
							dsnps+=snpcountchange;
							
							recsnps+=snphd;
						}
						
						// INCLUDE NEW ALLELE INTO BAC AND BAC INTO POP
						
						shared_ptr<Bac> RecBac = make_shared<Bac>(*bacs[recbac]);
						RecBac->genes[recgene]=RecGene;
						RecBac->emergeneration=Superpop::thisgeneration;
						// DO THE MAIN SEQUENCE RECOMBINATION AND REST OF METADATA ADJUSTING
						
						RecBac->recombine(recgene, rocus, donseq, snpcountchange, parand);
						
						// INCLUDE IN NEW STRAIN INTO POP
						
						bacs[recbac]=RecBac;
						
												 //  if (delhd> 0) {
												
													
													//   cout << "  AFTER AFTER AFTER AFTER AFTER AFTER AFTER AFTER AFTER AFTER  \n";

													//   BacRangeInfo(recbac, recbac, recgene);
													
													//   cout << "END END END END END END END END END END END END END END END END \n\n \n\n \n\n";
													// do{
													// 	cout << '\n' << "Press a key to continue...";
													// } while (cin.get() != '\n');

												 //  }
											  
						
						accrecs++; // increment accepted recombinations count
						Superpop::recs[(int)(hd*(donseq.length()))]++;
					
				
				} 
				else {
					Superpop::recs[-1]++;
				}	 	
			}
		}
		nextrec:
		;
	}		
	
	nrecombinations+=accrecs;
	
	
	return accrecs;
}


int Pop::recombineGene(int nrecs) { // Whole loci recombination
	
	//Variables for:
	int recbac=0; // index of recipient bacterium
	int donbac=0; //index of donor bacterium
	int recgene=0; // index of recombination gene
	int accrecs=0; // Amount of successful recombination tries.
	nrectries+=nrecs;
	
	vector<pair< int, string>>::iterator it;
	vector<pair< int, string>> difs;
	int inshd;
	
	vector<int>::iterator dit;
	vector<int> deldifs;
	
	int delhd=0;
	
	map< int, char>::iterator sit;
	map< int, char>::iterator sit2;
	
	int snphd=0;
	
	for (int i=0; i<nrecs; i++) {		
		donbac = bacs.size()*parand.getZeroone();
		recbac = bacs.size()*parand.getZeroone();
		recgene = ngenes*parand.getZeroone();
		
		if (bacs[donbac]->genes[recgene].get()==bacs[recbac]->genes[recgene].get()) { // Skip recombination actions if same allele
			nullrecombs++; 
			Superpop::recs[0]++;	
		} 
		else { 
			
			double tryrec=parand.getZeroone();
			double hd=0; 
				
			int totalsnps=bacs[donbac]->genes[recgene]->snps.size()+bacs[recbac]->genes[recgene]->snps.size();
			int totalindels=bacs[donbac]->genes[recgene]->dels.size()+bacs[recbac]->genes[recgene]->dels.size() + bacs[donbac]->genes[recgene]->inslens+bacs[recbac]->genes[recgene]->inslens;
			
			if (totalsnps==0 && totalindels==0) { // Skip recombination actions if original alleles
				nullrecombs++;
				Superpop::recs[0]++;
			}
			
			else {		// Recombination now with indel-extravaganza, hold your hat.
				
				
				// BEGIN Hamming Distance COMPUTING
				
				// hd from SNPS
				
				snphd=0;
				
				if (bacs[donbac]->genes[recgene]->snps.size()==0 || bacs[recbac]->genes[recgene]->snps.size()==0){
					snphd+=bacs[donbac]->genes[recgene]->snps.size() + bacs[recbac]->genes[recgene]->snps.size();
				}
				else {
					snphd=totalsnps;
					map< int, char>::iterator sit=bacs[donbac]->genes[recgene]->snps.begin();
					map< int, char>::iterator sit2=bacs[recbac]->genes[recgene]->snps.begin();
					
					while (sit != bacs[donbac]->genes[recgene]->snps.end() && sit2 !=bacs[recbac]->genes[recgene]->snps.end()){
						if (sit->first < sit2->first) {
							++sit;
						} else if (sit2->first < sit->first) {
							++sit2;
						} else { // equal keys
							if (sit->second==sit2->second) {
								snphd-=2;
							}
							if (sit != bacs[donbac]->genes[recgene]->snps.end()) {
								++sit;
							}
							if (sit2 !=bacs[recbac]->genes[recgene]->snps.end()){
								++sit2;
							}
						}
						
						if (sit == bacs[donbac]->genes[recgene]->snps.end()) {
							while (sit2 !=bacs[recbac]->genes[recgene]->snps.end()){
								if (sit->second==sit2->second) {
								snphd-=2;
								}
								++sit2;
							}
						} else if (sit2 ==bacs[recbac]->genes[recgene]->snps.end()){
							while (sit !=bacs[donbac]->genes[recgene]->snps.end()){
								if (sit->second==sit2->second) {
								snphd-=2;
								}
								++sit;
							}
						} 
						

					}
				}
				// END OF Hamming Distance COMPUTING
				
				hd+=snphd;  // Finally combine different hd counts
				
				
				// first hd from insertions
				if (totalindels >0 ) {
					
					if (tryrec > (double)(pow(10,-(hd/parand.genelen)*parand.recacc))) {
						Superpop::recs[-1]++;	
						goto nextrec;
					}	
					
					inshd=0;
					delhd=0;
					
					difs.resize(bacs[donbac]->genes[recgene]->ins.size()+bacs[recbac]->genes[recgene]->ins.size()); // Vector for differing insertions
					
					it=set_symmetric_difference(bacs[donbac]->genes[recgene]->ins.begin(),bacs[donbac]->genes[recgene]->ins.end(),bacs[recbac]->genes[recgene]->ins.begin(),bacs[recbac]->genes[recgene]->ins.end(), difs.begin()); // Get differing insertions, it points to last differing ins, i.e. last non null difs element.
					
					
					
					if (it!=difs.begin()) {  // If there is differences...
					difs.resize(it-difs.begin());			// drop null elements
					
					bool lastunique=1;
					
						for (int j=0; j<(int)difs.size()-1; j++ ) {
							
							if (difs[j].first==difs[j+1].first) { // if there is two inertions with same starting point (ins[n].first) compare.
								
								// Go through the string until either ends...
								for (int k= 0; k<(int)min(difs[j].second.length(),difs[j+1].second.length()); k++) {
									if (difs[j].second[k] != difs[j+1].second[k]) {
										inshd++;
									}
									
								}
								inshd+=fabs(difs[j].second.length()-difs[j+1].second.length()); // ...and add the diffence in length
								
								if (j+1==(int)difs.size()-1){
										lastunique=0;
								}
								
								j++;
							}
							
							else { // No insertion at same start point, add length of insertion to hd.
								inshd+=difs[j].second.length();
							}


							// If deemed preferable in the future, implementation could be added to reduce length of recombination in original sequence as more insertions are donated. 
							// if (difs[j].first <rocus+reclength-inshd) { 
							//		 break;
							// }
						}
						if (lastunique) {
							inshd+=difs[difs.size()-1].second.length();
						}
					}
					
					hd+=inshd;
					if (tryrec > (double)(pow(10,-(hd/parand.genelen)*parand.recacc))) {
						Superpop::recs[-1]++;	
						goto nextrec;
					}	
					// hd from deletion
					
					
					deldifs.resize(bacs[donbac]->genes[recgene]->dels.size()+bacs[recbac]->genes[recgene]->dels.size()); // Vector for differing deletions
					
					dit=set_symmetric_difference(bacs[donbac]->genes[recgene]->dels.begin(),bacs[donbac]->genes[recgene]->dels.end(),bacs[recbac]->genes[recgene]->dels.begin(),bacs[recbac]->genes[recgene]->dels.end(), deldifs.begin()); // Get differing insertions, it points to last differing ins, i.e. last non null difs element.
					
					delhd=distance(deldifs.begin(),dit);
					
					hd+=delhd;
				
				}
				
				hd=(double)(hd/(parand.genelen));  // Adjust with length
				
							
				if (hd==0) { // Skip recombination actions if original alleles
					nullrecombs++;
					Superpop::recs[0]++;
				}
				
				else if (tryrec<(double)(pow(10,-hd*parand.recacc))) { 

					if (parand.recdocthresh>0 && parand.recdocthresh <= hd*(parand.genelen)) {
						
						ostringstream filename;
						filename  << parand.fna;
						string fnadd=filename.str();
						
						ofstream recsfile; 
						recsfile.open ("outputs//RecEvents"+fnadd+".csv", ios::out | ios::app);
						recsfile << hd*(parand.genelen) <<";"<<(parand.genelen) << ";" << 0 << ";" << recgene  << ";" << bacs[donbac]->genes[recgene] << ";" << bacs[donbac] << ";" << bacs[recbac]<<";"<< popindex <<";"<< Superpop::thisgeneration <<endl ;
						recsfile.close();
					}
					
					
					
					int snpcountchange=bacs[donbac]->genes[recgene]->snps.size()-bacs[recbac]->genes[recgene]->snps.size();
					
					recsnps+=snphd;
					dsnps+=snpcountchange;
					
					// INCLUDE IN NEW ALLELE INTO BAC AND BAC INTO POP
					
					shared_ptr<Bac> RecBac = make_shared<Bac>(*bacs[recbac]);
					RecBac->genes[recgene]=bacs[donbac]->genes[recgene];

					if (parand.mlstdo>0) {
						RecBac->MLST=Superpop::getMLST(RecBac->genes);
					}
					
					if (parand.recsites==1) {
						RecBac->genes[recgene]->recsites.clear();
						RecBac->genes[recgene]->recsites[parand.genelen-1]=0;
					}
					
					RecBac->emergeneration=Superpop::thisgeneration;
					// DO THE MAIN SEQUENCE RECOMBINATION AND REST OF METADATA ADJUSTING
					
					// INCLUDE IN NEW STRAIN INTO POP
					
					bacs[recbac]=RecBac;
					
					accrecs++; // increment accepted recombinations count
					Superpop::recs[(int)(hd*parand.genelen)]++;

				}
				else {
					Superpop::recs[-1]++;
				}	
			} 
		}
		nextrec:
		;
	}		
	
	nrecombinations+=accrecs;
	
	
	return accrecs;
}



int Pop::insertions(int nins) { 
	int xbac=0;
	
	for (int i=0; i<nins; i++) {	 
		
		xbac = bacs.size()*parand.getZeroone();
		shared_ptr<Bac> xBac = make_shared<Bac>(*bacs[xbac]);
		
		xBac->insertion(parand);
		xBac->emergeneration=Superpop::thisgeneration;
		bacs[xbac]=xBac;
	}	
	return 0;
}

int Pop::deletions(int ndels) { 
	int xbac=0;
	
	for (int i=0; i<ndels; i++) {	
		
		xbac = bacs.size()*parand.getZeroone();
		shared_ptr<Bac> xBac = make_shared<Bac>(*bacs[xbac]);
		
		xBac->deletion(parand);
		xBac->emergeneration=Superpop::thisgeneration;
		bacs[xbac]=xBac;
	}	
	
	
	
	
	return 0;	
}


string Pop::migrate(Pop* dest, int migamount) { // Prompts migration from this population to dest.
	string migs="";
	int bactomigrate=0;
	map<int,int> temp;
	
	for (int i=0; i<migamount; i++) { 
		bactomigrate=bacs.size()*parand.getZeroone();
		int bactoreplace=parand.getZeroone()*(int)dest->bacs.size();	

		dest->bacs[bactoreplace]=bacs[bactomigrate];

		if (parand.migdo==1) {
			temp[bacs[bactomigrate]->MLST]++;
		}
	}
	
	bacsemigrated+=migamount;
	dest->bacsimmigrated+=migamount;
	
	if (parand.migdo==1) {
		for (auto const &it : temp) {
			ostringstream oss; 
			oss<< it.first <<"-" <<it.second << ",";
			migs+= oss.str();
		}
	}
	
	
	return migs;
}

string Pop::mepidemic() { 

	nmepidems++;
	int mepibac=bacs.size()*parand.getZeroone();
	int mepisize=parand.getPoisrand(parand.mepiSize);
	nmpedimbacs+=mepisize;
	int replacebac;
	string mics="";
	if (parand.micdo==1) {
		ostringstream oss; 
		oss << bacs[mepibac]->MLST <<";" << mepisize;
		mics= oss.str();
	}
	
	
	for (int i=0; i < mepisize; i++) {
		replacebac=bacs.size()*parand.getZeroone();

		bacs[replacebac]=bacs[mepibac];

	} 
	
	return mics;
}

void Pop::select() { 

	vector<shared_ptr<Bac>> temp ;	
	temp.reserve(bacs.size());

	int bacselect=0;
	for (int i=0; i < (int)bacs.size(); i++) {
		bacselect=bacs.size()*parand.getZeroone();
		temp.push_back(bacs[bacselect]);
	} 

	
	bacs=temp;  
	
		
}



void Pop::GeneWiseDiffs(const vector<int>  & randbacinds) {

    ostringstream diffs; 
    ostringstream mutdiffs; 

	
	int sharedlen=0;
	int alldifs=0;
	int mutationlydifs=0;
	
    int totaldiffs=0;
    int poptotaldiffs=0;

    int totalmutdiffs=0;
    int poptotalmutdiffs=0;
    double rm=0;
    
    map<pair<shared_ptr<Bac>,shared_ptr<Bac>>, pair<int,string>>  pairs; // Collection of different bacteria pairs (based on pointer addresses) with pair containing amount of pairs with the particular bacteria and differences etc. of them in a string.

    map<pair<shared_ptr<Bac>,shared_ptr<Bac>>, pair<int, string>>  mutpairs; 
    map<pair<shared_ptr<Bac>,shared_ptr<Bac>>, pair<int, int>>  muttotals; 
    
    forward_list<pair<shared_ptr<Bac>,shared_ptr<Bac>>> orderedpairs;
    
    for (int i =0 ; i< (int)randbacinds.size()-1; i++ ) { // Choose first bacterium (x) for pair
        for (int j = i+1; j< (int)randbacinds.size(); j++) { // Choose second bacterium (y) for pair
			
			pair<shared_ptr<Bac>,shared_ptr<Bac>> bacpair;
			
			if ( bacs[randbacinds[i]]< bacs[randbacinds[j]]) {
				bacpair= make_pair(bacs[randbacinds[i]], bacs[randbacinds[j]]);
			}
			else {
				bacpair = make_pair(bacs[randbacinds[j]], bacs[randbacinds[i]]);
			}
			
			
            auto it = begin(pairs);
            it = pairs.find(bacpair);
            
			if (it!=pairs.end()) {
				
				it->second.first=it->second.first+1;
				
			}
			
			else {
				ostringstream diffs; 

				ostringstream mutdiffs; 
				//~ cout << "Bacs" << i << "," << j << endl;
				
				orderedpairs.push_front(bacpair);
				
				if (bacs[randbacinds[i]].get() != bacs[randbacinds[j]].get()) { /* If the pair has different strains. .get() returns the address of the pointer. Some pointers may point to different addressess but actually contain same strain, but no two strains have same pointer. */
					
					totaldiffs=0;
					totalmutdiffs=0;
					for (int k =0; k<ngenes; k++) {
						
						if (bacs[randbacinds[i]]->genes[k].get() != bacs[randbacinds[j]]->genes[k].get()) { // Different gene pointer -> different allele (at least very probably)

							sharedlen=0;
							alldifs=0;
							mutationlydifs=0;
							
							GeneDistance(bacs[randbacinds[i]]->genes[k], bacs[randbacinds[j]]->genes[k], alldifs, sharedlen, mutationlydifs);
							
							diffs << ";"<< alldifs << ";" << sharedlen ;
							
							totaldiffs+=alldifs ;

							if (parand.mutonly==1) {
								mutdiffs << ";"<< mutationlydifs ;
								totalmutdiffs+=mutationlydifs ;
							}
							
							  //~ cout << "\n\n\n " ;
										//~ BacRangeInfo(randbacinds[i],randbacinds[i], k);
										//~ cout << "\n " ;
										//~ cout << "\n regular diffs" << Ijdif.size()+Jidif.size() << endl;
										//~ cout << "\n\n\n regular mutonly diffs" << mutationlydifs << endl;
										//~ cout << "\n ins diffs" << inshd << endl;
										//~ cout << "\n ins mutonly diffs" << mutationlyinsdifs<< endl;
										//~ cout << "\n " ;
										//~ cout << "\n i inslens " << bacs[randbacinds[i]]->genes[k]->inslens;
										//~ cout << "\n genelength+sharedinslen-pairdellens " << genelength+sharedinslen-pairdellens;
										//~ cout << "\n " ;
										//~ cout << "\n sharedinslen " << sharedinslen ;
										//~ cout << "\n " ;
										//~ cout << "\n pairdellens " << pairdellens ;
										//~ cout << "\n " ;
										//~ cout << "\n " ;
										//~ BacRangeInfo(randbacinds[j],randbacinds[j], k);
										//~ cout << "\n\n\n " ;

						}

						else { // same gene pointer -> same allele

							diffs << ";"<< 0  << ";" <<genelength+bacs[randbacinds[i]]->genes[k]->inslens-bacs[randbacinds[i]]->genes[k]->dels.size();
							if (parand.mutonly==1) {
								mutdiffs << ";"<< 0 ;
							}

						}
	 //END OF GENE COMPARISON                  
					}
				}
				else { // Same bac pointer -> same strain
					
					totaldiffs=0;

					if (parand.mutonly==1) {
						totalmutdiffs=0;
					}
					
					for (int k =0; k<ngenes; k++) {
						diffs << ";0;"<<genelength+bacs[randbacinds[i]]->genes[k]->inslens-bacs[randbacinds[i]]->genes[k]->dels.size() ;
						if (parand.mutonly==1) {
							mutdiffs << ";0" ;
						}

					}
					
				}

	//SUM WRITING
				diffs <<";" << totaldiffs  ;

				if (parand.mutonly==1) {
					muttotals[bacpair]={totaldiffs,totalmutdiffs};
				}
				
				pairs[bacpair]={1,diffs.str()};
				mutpairs[bacpair]={1,mutdiffs.str()};
			}
        }
    }
	
    ostringstream filename;
    filename << Superpop::thisgeneration << "-";

    filename << "Pop"<<popindex<<"-";
    

    filename   << parand.fna;
    string fnadd=filename.str();

    ofstream dfile; 
    dfile.open ("outputs//PairDistances"+fnadd+".csv");

    dfile<< "Bac_i_Bac_j;Amount_of_pairs";
    for (int i =0; i<ngenes; i++) {
        dfile << ";Gene_"<< i<< "_diffs"<< ";Shared_length_Gene_"<< i;
    }
    dfile << ";Total_diffs";
    dfile << endl;
    
	// for (auto & t : pairs) {
	// 	dfile<< t.first.first << "-" << t.first.second << ";" << t.second.first  << t.second.second << endl;
	// }

	for (auto & ot : orderedpairs) {
		auto t = pairs.find(ot);
		dfile<< t->first.first << "-" << t->first.second << ";" << t->second.first  << t->second.second << endl;
	}
	
    dfile.close();

    if (parand.mutonly==1) {
		ofstream mdfile; 
		mdfile.open ("outputs//MutationPairDistances"+fnadd+".csv");
		
		mdfile<< "Bac_i.Bac_j;Amount_of_pairs";
		for (int i =0; i<ngenes; i++) {
			mdfile << ";Gene_"<< i<< "_diffs";
		}
		mdfile << ";Total_diffs;Total_mutation_only_diffs;Population_diffs;Population_mutation_diffs;Population_R/M";
		mdfile << endl;
       
		for (auto & t : mutpairs) {
			 
			auto it2 = begin(pairs);
			it2 = pairs.find(t.first);
			
			mdfile<< t.first.first << "-" << t.first.second  << ";" << it2->second.first  << t.second.second;
			
			auto it = muttotals.find(t.first);
			 
			poptotaldiffs+= it2->second.first*it->second.first;
			poptotalmutdiffs+= it2->second.first*it->second.second;
			if (poptotalmutdiffs==0) {
				rm=0;
			}
			else {
				rm=(double)(poptotaldiffs-poptotalmutdiffs)/(poptotalmutdiffs);
			}
			mdfile <<";"<< it->second.first <<";" << it->second.second << ";" << poptotaldiffs << ";" << poptotalmutdiffs << ";" << rm<<endl;
		}
		mdfile.close();
    }


}


void Pop::InterpopGeneWiseDiffs(const vector<int>  & randbacinds,const vector<int>  & Brandbacinds, const Pop & popB) {
	
	int sharedlen=0;
	int alldifs=0;
	int mutationlydifs=0;
	
    int totaldiffs=0;
    int poptotaldiffs=0;

    int totalmutdiffs=0;
    int poptotalmutdiffs=0;
    double rm=0;
    
    map<pair<shared_ptr<Bac>,shared_ptr<Bac>>, pair<int,string>>  pairs; // Collection of different bacteria pairs (based on pointer addresses) with pair containing amount of pairs with the particular bacteria and differences etc. of them in a string.

    map<pair<shared_ptr<Bac>,shared_ptr<Bac>>, pair<int, string>>  mutpairs; 
    map<pair<shared_ptr<Bac>,shared_ptr<Bac>>, pair<int, int>>  muttotals; 

    forward_list<pair<shared_ptr<Bac>,shared_ptr<Bac>>> orderedpairs;
    
    for (int i =0 ; i< (int)randbacinds.size(); i++ ) { // Choose first bacterium (x) for pair
        for (int j = 0; j< (int)Brandbacinds.size(); j++) { // Choose second bacterium (y) for pair
			
			pair<shared_ptr<Bac>,shared_ptr<Bac>> bacpair;

			if ( bacs[randbacinds[i]]< popB.bacs[Brandbacinds[j]]) {
				bacpair= make_pair(bacs[randbacinds[i]], popB.bacs[Brandbacinds[j]]);
			}
			else {
				bacpair = make_pair(popB.bacs[Brandbacinds[j]], bacs[randbacinds[i]]);
			}
			
			
            auto it = begin(pairs);
            it = pairs.find(bacpair);
            
			if (it!=pairs.end()) {
				
				it->second.first=it->second.first+1;
				
			}
			
			else {
				ostringstream diffs; 

				ostringstream mutdiffs; 
				//~ cout << "Bacs" << i << "," << j << endl;
				orderedpairs.push_front(bacpair);
				// cout << bacpair.first << "-"<< bacpair.second << endl;
				
				if (bacs[randbacinds[i]].get() != popB.bacs[Brandbacinds[j]].get()) { /* If the pair has different strains. .get() returns the address of the pointer. Some pointers may point to different addressess but actually contain same strain, but no two strains have same pointer. */
					
					totaldiffs=0;
					totalmutdiffs=0;
					for (int k =0; k<ngenes; k++) {
						
						if (bacs[randbacinds[i]]->genes[k].get() != popB.bacs[Brandbacinds[j]]->genes[k].get()) { // Different gene pointer -> different allele (at least very probably)

							sharedlen=0;
							alldifs=0;
							mutationlydifs=0;
							
							GeneDistance(bacs[randbacinds[i]]->genes[k], popB.bacs[Brandbacinds[j]]->genes[k], alldifs, sharedlen, mutationlydifs);
							
							diffs << ";"<< alldifs << ";" << sharedlen ;
							
							totaldiffs+=alldifs ;

							if (parand.mutonly==1) {
								mutdiffs << ";"<< mutationlydifs ;
								totalmutdiffs+=mutationlydifs ;
							}
							
							  // cout << "\n\n\n " ;
									// 	BacRangeInfo(randbacinds[i],randbacinds[i], k);
									// 	cout << "\n " ;
									// 	popB.BacRangeInfo(Brandbacinds[j],Brandbacinds[j], k);
									// 	cout << "\n " ;
									// 	cout << "\n  diffs" << totaldiffs << endl;
									// 	cout << "\n  mutonly diffs" << mutationlydifs << endl;
									// 	cout << "\n\n\n " ;

						}

						else { // same gene pointer -> same allele

							diffs << ";"<< 0  << ";" <<genelength+bacs[randbacinds[i]]->genes[k]->inslens-bacs[randbacinds[i]]->genes[k]->dels.size();
							if (parand.mutonly==1) {
								mutdiffs << ";"<< 0 ;
							}

						}
	 //END OF GENE COMPARISON                  
					}
				}
				else { // Same bac pointer -> same strain
					
					totaldiffs=0;

					if (parand.mutonly==1) {
						totalmutdiffs=0;
					}
					
					for (int k =0; k<ngenes; k++) {
						diffs << ";0;"<<genelength+bacs[randbacinds[i]]->genes[k]->inslens-bacs[randbacinds[i]]->genes[k]->dels.size() ;
						if (parand.mutonly==1) {
							mutdiffs << ";0" ;
						}

					}
					
				}

	//SUM WRITING
				diffs <<";" << totaldiffs  ;

				if (parand.mutonly==1) {
					muttotals[bacpair]={totaldiffs,totalmutdiffs};
				}
				
				pairs[bacpair]={1,diffs.str()};
				mutpairs[bacpair]={1,mutdiffs.str()};
			}
        }
    }
	
    ostringstream filename;
    filename << Superpop::thisgeneration << "-";
    filename << "Pop"<<popindex<<"-"<< "Pop"<<popB.popindex<<"-";
    filename   << parand.fna;

    string fnadd=filename.str();

    ofstream dfile; 
    dfile.open ("outputs//InterpopPairDistances"+fnadd+".csv");

    dfile<< "Bac_i_Bac_j;Amount_of_pairs";
    for (int i =0; i<ngenes; i++) {
        dfile << ";Gene_"<< i<< "_diffs"<< ";Shared_length_Gene_"<< i;
    }
    dfile << ";Total_diffs";
    dfile << endl;
    
	for (auto & ot : orderedpairs) {
		auto t = pairs.find(ot);
		dfile<< t->first.first << "-" << t->first.second << ";" << t->second.first  << t->second.second << endl;
	}
	
    dfile.close();

    if (parand.mutonly==1) {
		ofstream mdfile; 
		mdfile.open ("outputs//InterpopMutationPairDistances"+fnadd+".csv");
		
		mdfile<< "Bac_i.Bac_j;Amount_of_pairs";
		for (int i =0; i<ngenes; i++) {
			mdfile << ";Gene_"<< i<< "_diffs";
		}
		mdfile << ";Total_diffs;Total_mutation_only_diffs;Population_diffs;Population_mutation_diffs;Population_R/M";
		mdfile << endl;
       
		for (auto & t : mutpairs) {
			 
			auto it2 = begin(pairs);
			it2 = pairs.find(t.first);
			
			mdfile<< t.first.first << "-" << t.first.second  << ";" << it2->second.first  << t.second.second;
			
			auto it = muttotals.find(t.first);
			 
			poptotaldiffs+= it2->second.first*it->second.first;
			poptotalmutdiffs+= it2->second.first*it->second.second;
			if (poptotalmutdiffs==0) {
				rm=0;
			}
			else {
				rm=(double)(poptotaldiffs-poptotalmutdiffs)/(poptotalmutdiffs);
			}
			mdfile <<";"<< it->second.first <<";" << it->second.second << ";" << poptotaldiffs << ";" << poptotalmutdiffs << ";" << rm<<endl;
		}
		mdfile.close();
    }


}


double Pop::Burner(const vector<int>  & randbacinds) {


	int nomutonlydifs=-1; 
	int sharedlen=0;
	int alldifs=0;
	
    int totaldiffs=0;
    unsigned long poptotaldiffs=0;

 
    map<pair<shared_ptr<Bac>,shared_ptr<Bac>>, pair<int, int>>  pairs; // Collection of different bacteria pairs (based on pointer addresses) with pair containing amount of pairs with the particular bacteria and total distance
	
    
    for (int i =0 ; i< (int)randbacinds.size()-1; i++ ) { // Choose first bacterium (x) for pair
        for (int j = i+1; j< (int)randbacinds.size(); j++) { // Choose second bacterium (y) for pair
			
			pair<shared_ptr<Bac>,shared_ptr<Bac>> bacpair;
			
			
			if ( bacs[randbacinds[i]]< bacs[randbacinds[j]]) {
				bacpair= make_pair(bacs[randbacinds[i]], bacs[randbacinds[j]]);
			}
			else {
				bacpair = make_pair(bacs[randbacinds[j]], bacs[randbacinds[i]]);
			}
			
			
            auto it = begin(pairs);
            it = pairs.find(bacpair);
            
			if (it!=pairs.end()) {
				
				it->second.first=it->second.first+1;
				
			}
			
			else {

				
				if (bacs[randbacinds[i]].get() != bacs[randbacinds[j]].get()) { /* If the pair has different strains. .get() returns the address of the pointer. Some pointers may point to different addressess but actually contain same strain, but no two strains have same pointer. */
					
					totaldiffs=0;
					
					for (int k =0; k<ngenes; k++) {
						
						if (bacs[randbacinds[i]]->genes[k].get() != bacs[randbacinds[j]]->genes[k].get()) { // Different gene pointer -> different allele (probably but not necessarily)

							sharedlen=0;
							alldifs=0;
							
							GeneDistance(bacs[randbacinds[i]]->genes[k], bacs[randbacinds[j]]->genes[k], alldifs, sharedlen, nomutonlydifs);
														
							totaldiffs+=alldifs ;

						}
                
					}
				}
				else { // Same bac pointer -> same strain
					totaldiffs=0;
				}

				pairs[bacpair]={1,totaldiffs};
				
			}
        }
    }
    
	int npairs=pairs.size();
	int nullpairs=0;
	
	for (auto & t : pairs) {
		if ( t.second.second >0 ) {
			poptotaldiffs+=  t.second.second ;
		}
		 else {
			nullpairs++;
		}
	}
	
	
	return (((double)poptotaldiffs/(npairs-nullpairs))/(double)(parand.genelen*parand.ngenes));
}



void Pop::GeneDistance(const shared_ptr<Gene> & igene, const shared_ptr<Gene> & jgene, int & alldifs, int & sharedlen, int & mutationlydifs) {  // Compares the given alleles at given gene of two bacs
	
	bool domutonly=parand.mutonly;
	if(mutationlydifs<0) {
		// To override parameter of mutation only difference counting
		domutonly=0;
		
	}

    int inshd=0;
    int sharedinslen=0;


    int mutationlyinsdifs=0;


    int maxdiffers = igene->snps.size() +jgene->snps.size();
    //MAIN SEQ SNPs
    vector<pair<int, char>> Ijdif(maxdiffers); // Get SNPs of bac_i not in bac_j
    vector<pair<int, char>>::iterator it;

    it=set_difference (igene->snps.begin(),igene->snps.end(), jgene->snps.begin(), jgene->snps.end(), Ijdif.begin());

    Ijdif.resize(it-Ijdif.begin());

    int pairdellens=igene->dels.size()+jgene->dels.size();

    if (igene->dels.size()>0) { // Erase SNPs aligned with deletion

        size_t delpos=igene->genstr.find("-");  // Decrease overlaps of deletion lengths
        while (delpos!=std::string::npos) {
            while (igene->genstr[delpos]=='-') {
                if (jgene->genstr[delpos]=='-') {
                    pairdellens--;
                }
                delpos++;
            }
            delpos=igene->genstr.find("-", delpos);
        }

        if (Ijdif.size() >0  ) {
            //~ cout << " del cleanup" << endl;

            //~ for (auto& t : Ijdif ) {
            //~ cout << t.second << ", " << t.first << endl;
            //~ }
            //~

            //~ cout << Ijdif.size() << endl;
            Ijdif.erase(remove_if(Ijdif.begin(), Ijdif.end(), [&](pair<int,char> & snpit ) {
                return jgene->genstr[snpit.first]=='-';
            }), Ijdif.end());



            //~ cout << Ijdif.size() << endl;
            //~ cout << " \n " << endl;
            //~ cout << igene->genstr  << endl;
            //~ cout << jgene->genstr  << endl;
            //~ cout << " \n " << endl;
            //~
            //~ for (auto& t : Ijdif ) {
            //~ cout << t.second << ", " << t.first << endl;
            //~ }
            //~
            //~ do
            //~ {
            //~ cout << '\n' << "Press a key to continue...";
            //~ } while (cin.get() != '\n');
            //~

            //~ cout << " \n\n" << endl;
        }
    }


    //~ cout << "The symmetric difference has " << (Ijdif.size()) << " elements:\n";
    //~ diffs << ";"<< (Ijdif.size()) ;
    //~ cout << igene->snps.size() << " " << jgene->snps.size() << endl;


    vector<pair<int, char>> Jidif(maxdiffers); // Get SNPs of bac_j not in bac_i

    it=set_difference (jgene->snps.begin(), jgene->snps.end(), igene->snps.begin(),igene->snps.end(), Jidif.begin());

    Jidif.resize(it-Jidif.begin());

    if (Jidif.size() >0  &&  jgene->dels.size()>0) { // Erase SNPs aligned with deletion

        //~ cout << " del cleanup" << endl;

        //~ for (auto& t : dif ) {
        //~ cout << t.second << ", " << t.first << endl;
        //~ }
        //~

        Jidif.erase(remove_if(Jidif.begin(), Jidif.end(), [&](pair<int,char> & snpit ) {
            return igene->genstr[snpit.first]=='-';
        }), Jidif.end());

    }

    if (domutonly==1) {
        mutationlydifs=0;

        //~ cout << "\nBac " << i << " gene " << k << endl;
        //~ for (auto& t : igene->recsites ) {
        // cout << t.second << ", " << t.first << endl;
        //~ }
        for (int o=0; o<(int)Ijdif.size(); o++) {
            // cout << dif[o].first << endl;

            map<int, int>::iterator itr;
            itr=igene->recsites.lower_bound(Ijdif[o].first);
            if (itr!=igene->recsites.end()) {
                if (itr->second > Ijdif[o].first) {
                    mutationlydifs++;
                }
                else {
                    //~ cout << " In recsite snp is " << itr->first << "  " << itr->second << endl;
                }
            }
            else {
                mutationlydifs++;
            }

        }
        //~ cout << "\nBac " << j << " gene " << k << endl;
        //~ for (auto& t : jgene->recsites ) {
        // cout << t.second << ", " << t.first << endl;
        //~ }

        for (int o=0; o<(int)Jidif.size(); o++) {
            // cout << Jidif[o].first << endl;

            map<int, int>::iterator itr;
            itr=jgene->recsites.lower_bound(Jidif[o].first);
            if (itr!=jgene->recsites.end()) {
                if (itr->second > Jidif[o].first) {
                    mutationlydifs++;
                }
                //~ else {
                //~ cout << " In recsite snp is " << itr->first << "  " << itr->second << endl;
                //~ }
            }
            else {
                mutationlydifs++;
            }

        }


    }

    sharedinslen=0;

    // If there is insertions, look for SNPs in the ones that align completely (Same start point and legnth).
    // This means that two different insertions can be counted to distance as well as SNP differnces.
    if (igene->ins.size()>0 && jgene->ins.size()>0) {
        // INSERT SNPs
        vector<pair< int, string>>::iterator itins;

        vector<pair< int, string>> Ijinsdifs(igene->ins.size()+jgene->ins.size()); // Vector for differing insertions

        vector<int> Isnpins;
        itins=set_difference(igene->ins.begin(), igene->ins.end(), jgene->ins.begin(), jgene->ins.end(), Ijinsdifs.begin()); // Get differing insertions, itins points to last differing ins,

        Ijinsdifs.resize(itins-Ijinsdifs.begin());

        vector<pair< int, string>> Jiinsdifs(igene->ins.size()+jgene->ins.size()); // Vector for differing insertions

        vector<int> Jsnpins;
        itins=set_difference(jgene->ins.begin(), jgene->ins.end(), igene->ins.begin(), igene->ins.end(), Jiinsdifs.begin()); // Get differing insertions, itins points to last differing ins,

        Jiinsdifs.resize(itins-Jiinsdifs.begin());
        inshd=0;

        sharedinslen=igene->inslens;

        if (Ijinsdifs.size()>0) {

            for (int d=0; d<(int) Ijinsdifs.size(); d++ ) {
                bool sharedins=0;
                for (int e=0; e<(int) Jiinsdifs.size(); e++ ) {

                    if ( Ijinsdifs[d].first==Jiinsdifs[e].first && Ijinsdifs[d].second.length()==Jiinsdifs[e].second.length() ) { // if there is two insertions with same starting point (insdifs[n].first) and length, compare.
                        sharedins=1;

                        for (int k= 0; k<(int)Jiinsdifs[e].second.length(); k++) {
                            if ( Ijinsdifs[d].second[k] !=  Jiinsdifs[e].second[k]) {
                                inshd++;
                                if (isupper(Jiinsdifs[e].second[k])) {
                                    Jsnpins.push_back(Jiinsdifs[e].first);
                                }
                                else if (isupper(Ijinsdifs[d].second[k])) {
                                    Isnpins.push_back(Ijinsdifs[d].first);
                                }
                            }
                        }
                    }
                }
                if (sharedins==0) {
                    sharedinslen-=Ijinsdifs[d].second.length();
                }

            }

            if (domutonly==1 && Jsnpins.size()+Isnpins.size()>0) {
                // insSNPs() return count of "real" SNPs not in recsite.
                mutationlyinsdifs=0;

                for (int o=0; o<(int)Isnpins.size(); o++) {
                    // cout << insdif[o].first << endl;
                    mutationlyinsdifs+=igene->recsiteSNPs(Isnpins[o]);

                }
                //~ cout << "^ " << randbacinds[i] << "      v" << randbacinds[j] << endl;
                for (int o=0; o<(int)Jsnpins.size(); o++) {
                    // cout << Jiinsdif[o].first << endl;

                    mutationlyinsdifs+=jgene->recsiteSNPs(Jsnpins[o]);

                }


            }
        }
        /*INSERTS END*/
    }


    alldifs=Ijdif.size()+Jidif.size()+inshd;
    mutationlydifs+=mutationlyinsdifs;
    sharedlen= genelength+sharedinslen-pairdellens;



}




string Pop::summary(int generation, int popindex) { 

	ostringstream sumstr; 
	sumstr.unsetf ( ios::skipws );
	
	sumstr << "                   Population: " << popindex  << " in generation: " << generation << endl;
	sumstr << "\nSize: " << bacs.size() << endl;
	
	sumstr << "\nNumber of mutations " << tmuts << endl;
	sumstr << "Number of mutations on previous SNPs " << mutsonsnp << endl;
	sumstr << "Number of mutations to ancestral from SNP " << mutstoanc << endl;
	sumstr << "\nNumber of recombinations " << nrecombinations+nullrecombs << endl;
	sumstr << "Number of null recombinations " << nullrecombs << ", of all recombinations "<< (double)nullrecombs/(nrecombinations+nullrecombs) << endl;
	sumstr << "Number of recombination tries " << nrectries << endl;
	
	sumstr << "\nNumber of snps introduced by mutations during the whole simulation " << snpsfromuts << endl;
	sumstr << "Number of snps added by recombinations during the whole simulation " << recsnps << endl;
	sumstr << "Ratio of the two above, #(snps_rec)/#(snps_muts) " << (double)recsnps/tmuts << endl;
	
	sumstr << "\nBacteria immigrated: " << bacsimmigrated << endl;
	sumstr << "Bacteria emigrated: " << bacsemigrated << endl;
	
	sumstr << "\nNumber of microepidemics: " << nmepidems << endl;
	sumstr << "\nNumber of bacteria from microepidemics: " << nmpedimbacs << endl;
	
	if (parand.ngenerations >generation) {
		sumstr << "\n  \n  " << endl;
	}
	summarystr += sumstr.str();
	
	return sumstr.str();
}


void Pop::PopInfo() {
	
	for (int i =0 ; i< (int)bacs.size(); i++ ) { 

				
				cout << " \n \n ----------- Bacs " <<  i << "   ptr-address: " << bacs[i].get()   << "   snpn " << bacs[i]->snpn <<  " -----------" << endl;

				cout << "Mutations: " << bacs[i]->nmutations << endl;
				cout << "Recombinations: " << bacs[i]->nrecombinations << endl;
				
				cout << " \nRecsites : " <<endl;
				for (int m=0; m<ngenes ; m++) {
					cout << " recsites size " << bacs[i]->genes[m]->recsites.size()<< endl;
					if (bacs[i]->genes[m]->recsites.size() > 0) {
						cout << m <<": ";
					
						for ( auto &k : bacs[i]->genes[m]->recsites ) {
							cout << k.second << ", " << k.first << " ";
						}
						cout <<"\n";
					}
				}
				
				cout <<"\n";
				
				cout << "snps, Bac" << i << endl;

				for (int m=0; m<ngenes ; m++) {
					cout << m<< ": " << "ptr-address: " << bacs[i]->genes[m].get() << " \n";
					for ( auto &k : bacs[i]->genes[m]->snps ) {
					cout << k.first << " " << k.second << ", ";
					}

					cout <<"\n\n";
				}
				
				cout <<"\n\n";

				bacs[i]->printGenes();

				
				
	}
		
	cout << " \n \n " << endl;
	
}

void Pop::BacRangeInfo(int lowbound, int upbound) {
	cout << "*************************************************************" << endl;
	upbound++;
	if(lowbound>upbound) {
		cout << "Lower bound given for BacRangeInfo() too big)" << endl;
	return;
	}
	
	if(upbound>(int)bacs.size()) {
		cout << "Upper bound given for BacRangeInfo() too big)" << endl;
	return;
	}
	
	
	for (int i =lowbound ; i< upbound; i++ ) { 

				cout << " \n \n ----------- Bacs " <<  i << "   ptr-address: " << bacs[i].get()   << "   snpn " << bacs[i]->snpn <<  " -----------" << endl;

				cout << "Mutations: " << bacs[i]->nmutations << endl;
				cout << "Recombinations: " << bacs[i]->nrecombinations << endl;
				
				cout << " \nRecsites : " <<endl;
				for (int m=0; m<ngenes ; m++) {
					cout << " recsites size " << bacs[i]->genes[m]->recsites.size()<< endl;
					if (bacs[i]->genes[m]->recsites.size() > 0) {
						cout << m <<": ";
					
						for ( auto &k : bacs[i]->genes[m]->recsites ) {
							cout << k.second << ", " << k.first << " ";
						}
						cout <<"\n";
					}
				}
				
				cout <<"\n";
				
				cout << "snps, Bac" << i << endl;

				for (int m=0; m<ngenes ; m++) {
					cout << m<< ": " << "ptr-address: " << bacs[i]->genes[m].get() << " \n";
					for ( auto &k : bacs[i]->genes[m]->snps ) {
					cout << k.first << " " << k.second << ", ";
					}

					cout <<"\n\n";
				}
				
				cout <<"\n\n";
				
	}
		
	cout << " \n \n " << endl;
	
	cout << "*************************************************************\n\n" << endl;
}

void Pop::BacRangeInfo(int lowbound, int upbound, int m) {
	cout << "*************************************************************" << endl;
	upbound++;
	if(lowbound>upbound) {
		cout << "Lower bound given for BacRangeInfo() too big)" << endl;
	return;
	}
	
	if(upbound>(int)bacs.size()) {
		cout << "Upper bound given for BacRangeInfo() too big)" << endl;
	return;
	}
	
	
	for (int i =lowbound ; i< upbound; i++ ) {

				cout << " \n ----------- Bacs " <<  i << "   ptr-address: " << bacs[i].get()   << "   snpn " << bacs[i]->snpn <<  " -----------" << endl;
				cout << "Mutations: " << bacs[i]->nmutations << endl;
				cout << "Recombinations: " << bacs[i]->nrecombinations << endl;
				
				
				cout << "Recsites size " << bacs[i]->genes[m]->recsites.size()<< endl;
				if (bacs[i]->genes[m]->recsites.size() > 0) {
				
					for ( auto &k : bacs[i]->genes[m]->recsites ) {
						cout << k.second << ", " << k.first << " ";
					}
					cout <<"\n";
				}
				 
				cout << "Insertions size " << bacs[i]->genes[m]->ins.size()<< endl;
				if (bacs[i]->genes[m]->ins.size() > 0) {
					
					for ( auto &k : bacs[i]->genes[m]->ins ) {
							cout << k.first << "  " << k.second.length() << " " << k.second << " ";
						
					}
					cout <<"\n";
				}
				
				 
				cout << "Ins-recsites size " << bacs[i]->genes[m]->insrecs.size()<< endl;
				if (bacs[i]->genes[m]->insrecs.size() > 0) {
				
					for ( auto &k : bacs[i]->genes[m]->insrecs ) {
						
						for ( auto &l : k.second ) {
							cout << k.first << "  " << l.second << "-" << l.first<< " ";
						}
					}
					cout <<"\n";
				}
				
				cout << "Deletions size " << bacs[i]->genes[m]->dels.size()<< endl;
				if (bacs[i]->genes[m]->dels.size() > 0) {
				
					for ( auto &k : bacs[i]->genes[m]->dels ) {
						
						cout << k <<  " ";
						
					}
					cout <<"\n";
				}

				cout << "Del lens " << bacs[i]->genes[m]->dels.size()<< endl;
				
				
				cout <<"\n";
				
				cout <<"\nGene " <<  m<< ": " << "ptr-address: " << bacs[i]->genes[m].get() << " \n";
				for ( auto &k : bacs[i]->genes[m]->snps ) {
				cout << k.first << " " << k.second << ", ";
				}
				cout << "\n"<< bacs[i]->genes[m]->getGeneText(1);
				
				
				
				cout <<"\n";
				
	}
		
	cout << " \n " << endl;
	
	cout << "*************************************************************\n" << endl;
}

int Pop::Diversity(int generation, int popindex) {
	
	ostringstream filename;
	
	filename << "Pop"<<popindex<<"-";
	
	
	filename   << parand.fna;
	string fnadd=filename.str();
    
	ofstream divs;
	divs.open ("outputs//DiversityDistribution"+fnadd+".csv", ios::out | ios::app);
	
	
	map<shared_ptr<Bac>, int> types;
	types.clear();

	for (int i =0 ; i < nbacs; i++){
		
		types[bacs[i]]++;
	}
	for (auto const &iter : types) {
		
		divs << generation << ";" << iter.first << ";"  << iter.first->emergeneration << ";" << iter.first->snpn << ";" <<  iter.second  << ";" << (double)iter.second/nbacs << endl;
		
	}
	
	divs.close();

	
	return types.size();	
}


