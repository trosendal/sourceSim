#include "Gene.h"
#include "Para.h"
#include "Superpop.h"
#include <random>
#include <string>
#include <sstream>
#include <algorithm>


using namespace std;

/* "Gene" class represents a single loci of given base pair length. It is the smallest level of this simulator. The term "gene" is used for historical reasons, and thus might be used interchangeably with loci in the code files or readme.

Gene creates and stores the sequence of loci in a string and variables for various uses. It also handles most of the loci-level events.

This class only uses one other class, which is the Para. Here the Para class is used by Para object named "pa". Para objects hold all the input parameters and handle random number generation.

The bases are always listed in same order in this simulator, for example in vectors. The order is:
0. A  
1. T  
2. G  
3. C
In retrospect this isn't best order. However, if it helps remembering the order, it follows from starting with "A" and having G and C in the end together as in "GC-ratio".


*/

double Gene::mutmatrix[12];

Gene::Gene(int genelen, Para pa, int ng)  {

	genei=ng;	
	genstr="";

	genelength=genelen; 
	inslens=0;

	nmutations=0;
	nrecombinations=0;

	allele=1;
	mlstdo=pa.mlstdo;
	recsitesdo= pa.recsites;

	updategeneration=0;

	aa=0; 
	cc=0;		
	gg=0;	
	tt=0;	

	for (int i=0; i<genelength; i++) {  
	
		
		double baserand=pa.getZeroone();  
		vector<double> bases;  
		bases.reserve(3);   
		bases=pa.getBaseprops();  
		
		if (baserand<bases[0]) {  
			genstr+="a";  
			aa++; 
		}
		else if (baserand <(bases[1]) && baserand>bases[0]) {  
			genstr+="t";  
			tt++;
		}  
		else if (baserand <(bases[2]) && baserand>(bases[1])) {  
			genstr+="g";
			gg++;
		}
		else {    
			genstr+="c";
			cc++;
		}
	}
	
	 
	mutmatrix[0]=pa.atot; mutmatrix[1]=pa.atog+pa.atot; mutmatrix[2]=pa.atoc+pa.atog+pa.atot;
	mutmatrix[3]=pa.ttoa; mutmatrix[4]=pa.ttog+pa.ttoa; mutmatrix[5]=pa.ttoc+pa.ttog+pa.ttoa;
	mutmatrix[6]=pa.gtoa; mutmatrix[7]=pa.gtot+pa.gtoa; mutmatrix[8]=pa.gtoc+pa.gtot+pa.gtoa;
	mutmatrix[9]=pa.ctoa; mutmatrix[10]=pa.ctot+pa.ctoa; mutmatrix[11]=pa.ctog+pa.ctot+pa.ctoa;
	
	
}


string Gene::getGeneText(bool inserts) { 	 
		
	if (inserts) {
		string tempstr=genstr;
		std::map<int, string>::reverse_iterator rit;
		for (rit=ins.rbegin(); rit!=ins.rend(); ++rit) {
			string tempins="("+rit->second+")";
			tempstr.insert(rit->first, tempins);
		}
		
		return tempstr;
	}
	else {
		return genstr;
	}
}


void Gene::mutate(Para &pa, int& snpd)  
{
	
	char oldbase='-';
	int position;
	position = (genelength+inslens)*pa.getZeroone();  
	
	if (position <genelength) {  
		
		oldbase = genstr[position];
		while (oldbase=='-') {
			position = genelength*pa.getZeroone();  
			oldbase = genstr[position];
		}
		
		double baserand=pa.getZeroone();

		char newbase = mutationBase(oldbase, baserand);
		

		//~ cout << "\nGene   snps \n" ;
		//~ for ( auto &k :snps ) {
		//~ 	cout << k.first  << " " << k.second << endl; 
		//~ }  
							
		nmutations++; 


		if (isupper(oldbase)){  // If mutation happened on a position that's already mutated....
				
			if (tolower(newbase)==Para::origenes[genei][position]) { // If that position mutated back to initial base...
				std::map< int, char>::iterator it = snps.find(position);
				
				if (it != snps.end()){
					snps.erase(it);
				}
				
				newbase=tolower(newbase);
				snpd=-1;
			}
					
			else { 									// or to new derived snp
				std::map< int, char>::iterator it = snps.find(position);
				if (it != snps.end()){
					it->second = newbase;
				}
				
				snpd=0;
			}	
		} 
		else {
			snps[position]=newbase;
			snpd=1;
		}
		
		genstr[position]=newbase;
		


		if (recsitesdo==1) {
			mutationRecsitesplit(position);
		}
		if (mlstdo>0) {
			allele=Superpop::getAllele(getGeneText(1), genei);
		}

	}
	
	
	else { // Mutation on insert-site
			
		map< int, string>::iterator it = ins.begin();
		int outofins=genelength;
		
		while (it!=ins.end()) {
			outofins=outofins+it->second.length();
			if (position < outofins) {
				position=fabs(position-outofins)-1;
				break;
			}
			else {
				it++;
			} 
		}
		
		
		char oldbase=it->second[position];
		
		while (oldbase=='-') {
			position = it->second.length()*pa.getZeroone();  
			oldbase = it->second[position];
		}
			
		double baserand=pa.getZeroone();  
			
		char newbase = mutationBase(oldbase, baserand);
		
							
		nmutations++;  

		if (isupper(oldbase)){   
			snpd=0;	
		} 
		else {
			snpd=1;
		}

		it->second[position]=newbase;
		

		
		if (recsitesdo==1) {
			mutationInsRecsitesplit(it->first,position);
		}
		if (mlstdo>0) {
			allele=Superpop::getAllele(getGeneText(1), genei);
		}
		
	
	}
}


char Gene::mutationBase(char oldbase, double baserand) {
	
		char newbase='W';
		
		switch ( oldbase ) {  
		case 'A':  
		case 'a':

			if (baserand<mutmatrix[0]) {  
				newbase='T';
			} 
			else if (baserand <(mutmatrix[1]) && baserand>mutmatrix[0]) { 
				newbase='G';
			} 
			else if (baserand <(mutmatrix[2]) && baserand>(mutmatrix[1])) { 
				newbase='C';
			} 
			break;
			
		case 't':
		case 'T':
			if (baserand<mutmatrix[3]) { 
				newbase='A';
			} 
			else if (baserand <(mutmatrix[4]) && baserand>mutmatrix[3]) { 
				newbase='G';
			} 
			else if (baserand <(mutmatrix[5]) && baserand>(mutmatrix[4])) { 
				newbase='C';
			} 
			break;
		case 'g':
		case 'G':
			if (baserand<mutmatrix[6]) { 
				newbase='A';
			} 
			else if (baserand <(mutmatrix[7]) && baserand>mutmatrix[6]) { 
				newbase='T';
			} 
			else if(baserand <(mutmatrix[8]) && baserand>(mutmatrix[7])) { 
				newbase='C';
			} 
			break;
		case 'c':
		case 'C':
			if (baserand<mutmatrix[9]) { 
				newbase='A';
			} 
			else if (baserand <(mutmatrix[10]) && baserand>mutmatrix[9]) { 
				newbase='T';
			} 
			else if (baserand <(mutmatrix[11]) && baserand>(mutmatrix[10])) { 
				newbase='G';
			} 
			break;	
		default:
			cout << "Base to be mutated unknown: " << oldbase << endl;
			break;
		}
	
	return newbase;
}


void Gene::recombine (int rocus, string donatedseq) {  
	
	int reclen=donatedseq.length();
	genstr.replace(rocus, reclen, donatedseq); 
	
	nrecombinations++; 
	 
	inslens=0; 

	for (auto & g: ins ) {
		inslens+=g.second.length();
			if (rocus<g.first && g.first<rocus+reclen) {
				int tempint=g.first;
				insrecs[tempint]=vector<pair<int,int>>{{g.second.length()-1 , 0}};
			}
	}			

	if (mlstdo>0) {
		allele=Superpop::getAllele(getGeneText(1), genei);
	}

	if (recsitesdo==1) {
		if (nrecombinations==0 ) {
			recsites[rocus+reclen]=rocus;
		}
		else {
			addRecsite(rocus+reclen, rocus);
		}
	}

}


void Gene::addRecsite (int endpoint, int startpoint) { /* Add new recsite entry from startpoint to endpoint,  with merging to others if applicable. 
This is done by finding smallest previous end points that are still bigger than endpoint and startpoint of entry recsite. 
Overlaps with entry site are merged so that greatest end point and smallest startpoint in the overlaps are preserved. 

Utilises maps where recsites[endpoint]=startpoint and recsiteX.first=endpoint, recsiteX.second=startpoint. 

Endpoint is used as key because of lower_bounds possible return of .end() */
	
	std::map<int,int>::iterator itlow,ithigh;
	
	itlow=recsites.lower_bound(startpoint);  // Get from previous recsites smallest endsite that is higher or equal to entry start.
			
			if (itlow==recsites.end()) { // No previous endpoints greater than startpoint, add new site.
				recsites[endpoint]=startpoint; // Entry recsite is after all previous recsites
				return;
			}
			
			ithigh=recsites.lower_bound(endpoint);  // Get from previous recsites the smallest end that's higher or equal to entry end.
				
				
			if (ithigh==recsites.end()) { // Entry startpoint smaller than some previous, entry endpoint is higher that previous ones.
					
				if (itlow->second < startpoint ) { // Entry recsite partially overlaps previous one(s). Merging from entry endpoint to previous startpoint.
						
						int temp=itlow->second; // Replace previous recsite with combination of the entry and the previous one.
						recsites.erase(itlow, recsites.end());
						recsites[endpoint]=temp;
						
				}
					
				else { // Entry recsite completely overlaps found one(s). Merging from entry endpoint to found startpoint.
						
						recsites.erase(itlow,recsites.end());
						recsites[endpoint]=startpoint;
				}
			}
			
			else { // Endpoint lower than previous endpoint
				
				if (ithigh->second < startpoint ) { // Entry site completely within previous one. Do nothing
				}
				else if ( ithigh == itlow) { // Entry startpoint between found start point and end point of the possible recsites before the found one.
					if (ithigh->second < endpoint ) { // Entry recsite partially overlaps previous one(s). Merging from previous endpoint to entry startpoint.
						ithigh->second=startpoint;
					}
					else {					
					// Entry site between previous ones. Add new recsite 
					recsites[endpoint]=startpoint;
					}
				}
				else { // Entry start point smaller than some previous endpoint, which is smaller than entry endpoint.
					if (ithigh->second < endpoint) { // Entry recsite overlaps on the end and startside. Replace them with merge of all.
						int temp = itlow->second;
						recsites.erase(itlow, ithigh);
						ithigh->second= temp;
					}
					else { // Only start side overlaps. Replace with merge
						int temp=itlow->second;
						recsites.erase(itlow, ithigh);
						recsites[endpoint]=temp;
					}
						
				}
			}

}


void Gene::mutationRecsitesplit(int position) {
 			
	if ((int)recsites.size()>0) {

		std::map<int,int>::iterator itlow;
		itlow=recsites.lower_bound(position);
		if (itlow!=recsites.end()) {
		
			if (itlow->second <= position ) {
				if (position==genelength) {
					recsites[position-1]=itlow->second;
					recsites.erase(itlow);
				}
				else if (position==0) {
					itlow->second=1;
				}
				else {
					recsites[position-1]=itlow->second;
					itlow->second=position+1;
				}
			}
		}
	
	}
}


void Gene::mutationInsRecsitesplit(int inskey, int position) {


    map<int, vector<pair<int,int>>>::iterator inss = insrecs.find(inskey);
    if (inss!=insrecs.end()) {
        vector<pair<int,int>>::iterator rsite=inss->second.begin();
        while (rsite!=inss->second.end()) {
            if (rsite->first-rsite->second<2) {
                inss->second.erase(rsite);
            }
            else if (rsite->second < position  &&  position < rsite->first) {
                int temp=rsite->second;
                rsite->second=position+1;
                inss->second.push_back(make_pair(position-1,temp));
                break;
            }
            else if (position == rsite->first) {
                rsite->first=position-1;
                break;
            }
            else if (position == rsite->second) {
                rsite->second=position+1;
                break;
            }
            else {
                rsite++;
            }
        }
    }
}

		
int Gene::recsiteSNPs(int inskey) {  // Get amount of snps in recsite. Input must be starting point of insertion.

    int mutsnps=0;

    map<int, vector<pair<int,int>>>::iterator recsites = insrecs.find(inskey);
    if (recsites!=insrecs.end()) {

        string insertion=ins.find(inskey)->second;

        for (int i=0; i<(int)insertion.length(); i++ ) {

            if (isupper(insertion[i])) {
				mutsnps++;
                vector<pair<int,int>>::iterator rsite=recsites->second.begin();
                while (rsite!=recsites->second.end()) {
                    if (rsite->second<=i && i<=rsite->first) {
                       mutsnps--;
                    }
                    rsite++;

                }
            }
        }
    }

    return mutsnps;
}


void Gene::insertion(Para &pa) {
	
	
	int len=0;
	
	while (len <1 || len > pa.indelmax*genelength) {
			len=pa.getZipf(0,1.0);
	}
	Superpop::indels[len]++;
	
	inslens+=len;
	
	vector<double> bases;  
	bases.reserve(3);   
	bases=pa.getBaseprops();  
	
	string newbases="";
	
	while (len--) {
		
		double baserand=pa.getZeroone();
		if (baserand<bases[0]) {  
			newbases+="a";  
		}
		else if (baserand <(bases[1]) && baserand>bases[0]) {  
			newbases+="t";  
		}  
		else if (baserand <(bases[2]) && baserand>(bases[1])) {  
			newbases+="g";
		}
		else {    
			newbases+="c";
		}
	}	
	
	int position=genelength*pa.getZeroone();
	
	while (genstr[position]=='-') {
		position = genelength*pa.getZeroone();  
	}
	
	pair<int,string> temppair (position, newbases);
	pair<map<int,string>::iterator,bool> inserttest=ins.insert(temppair);
	
	
	while(inserttest.second==false) {
		position=genelength*pa.getZeroone();
		while (genstr[position]=='-') {
			position = genelength*pa.getZeroone();  
		}
		temppair.first=position;
		inserttest=ins.insert(temppair);
	}

	if (mlstdo>0) {
		allele=Superpop::getAllele(getGeneText(1), genei);
	}
	
}


void Gene::deletion(int rocus, int len) {

	std::set<int>::iterator it=dels.begin();
	
	for ( int i=rocus; i<rocus+len; i++ ) {
				
		it = dels.insert(it,i);
			switch ( genstr[i] ) {  
					case 'A': 
					case 'a':
						aa--;
						break;
				case 'T': 
					case 't':
						tt--;  
						break;
				case 'G': 
					case 'g':
						gg--;  
						break;
				case 'C': 
					case 'c':
						cc--;  
						break;
				case '-':
						break;
				default:
					cout << "Nucleotide count error!" << endl;
				
			}
	}
	
	genstr.replace(rocus,len,len,'-');  
	
	
	if (snps.size()>0) {	

		std::map<int,char>::iterator itlow,itup;		
		itup=snps.upper_bound(rocus+len-1); 
		itlow=snps.lower_bound(rocus); 
		
		snps.erase(itlow,itup);       

	}	
		
	if (ins.size()>0) {	
		 
		std::map<int,string>::iterator inslow,insup;		
		insup=ins.upper_bound(rocus+len-1); 
		inslow=ins.lower_bound(rocus); 
		for (auto it = inslow; it != insup; ++it) { 
			inslens-=it->second.length();
		}
		ins.erase(inslow,insup);       
	
	}		
					
	if (mlstdo>0) {
		allele=Superpop::getAllele(getGeneText(1), genei);
	}
}


void Gene::updateBasecounts() {
		
	string oristr=Para::origenes[genei];
		for ( auto &k : snps ) {
				char snpbase=k.second;
				int snpposition=k.first;
				switch ( snpbase ) {  
					case 'A': 
					case 'a':
						aa++;  
						if ( oristr[snpposition]=='t') {  
							tt--;  
						} 
						else if (oristr[snpposition]=='g') { 
							gg--;
						} 
						else if (oristr[snpposition]=='c') { 
							cc--;
						} 
						break;
				case 'T': 
					case 't':
						tt++;  
	
						if ( oristr[snpposition]=='a') {  
							aa--;  
						} 
						else if (oristr[snpposition]=='g') { 
							gg--;
						} 
						else if (oristr[snpposition]=='c') { 
							cc--;
						} 
						break;
				case 'G': 
					case 'g':
						gg++;  
	
						if ( oristr[snpposition]=='a') {  
							aa--;  
						} 
						else if (oristr[snpposition]=='t') { 
							tt--;
						} 
						else if (oristr[snpposition]=='c') { 
							cc--;
						} 
						break;
				case 'C': 
					case 'c':
						cc++;  
	
						if ( oristr[snpposition]=='a') {  
							aa--;  
						} 
						else if (oristr[snpposition]=='t') { 
							tt--;
						} 
						else if (oristr[snpposition]=='g') { 
							gg--;
						} 
						break;
				default:
					cout << "Nucleotide count error!" << endl;
				
			}
		}
	
		for ( auto &k : ins ) {
				string bases=k.second;
				for (int i=0; i<(int)bases.length(); i++){
					
					switch ( bases[i] ) {  
						case 'A': 
						case 'a':
							aa++;  
		
							break;
					case 'T': 
						case 't':
							tt++;  
		
							break;
					case 'G': 
						case 'g':
							gg++;  
		
							
							break;
					case 'C': 
						case 'c':
							cc++;  
		
							
							break;
					default:
						cout << "Nucleotide count error!" << endl;
				}
			}
		}
	

}




void Gene::updbasecountbySNP(int snpposition, char snpbase) {
		
	char oribase=Para::origenes[genei][snpposition];
	
	switch ( snpbase ) {  
		case 'A': 
		case 'a':
			aa--;  
	
			if ( oribase=='t') {  
				tt++;  
			} 
			else if (oribase=='g') { 
				gg++;
			} 
			else if (oribase=='c') { 
				cc++;
			} 
			break;
	case 'T': 
		case 't':
			tt--;  
	
			if ( oribase=='a') {  
				aa++;  
			} 
			else if (oribase=='g') { 
				gg++;
			} 
			else if (oribase=='c') { 
				cc++;
			} 
			break;
	case 'G': 
		case 'g':
			gg--;  
	
			if ( oribase=='a') {  
				aa++;  
			} 
			else if (oribase=='t') { 
				tt++;
			} 
			else if (oribase=='c') { 
				cc++;
			} 
			break;
	case 'C': 
		case 'c':
			cc--;  
	
			if ( oribase=='a') {  
				aa++;  
			} 
			else if (oribase=='t') { 
				tt++;
			} 
			else if (oribase=='g') { 
				gg++;
			} 
			break;
	default:
		cout << "Nucleotide count error!" << endl;
	
	}
	
	

}

