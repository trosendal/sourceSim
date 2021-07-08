#include "Bac.h"
#include "Superpop.h"
#include "Para.h"
#include <random>
#include <string>
#include <sstream>
#include <memory>
#include <set>


using namespace std;
Bac::Bac(int ngene, int genelen, Para pa)  {   
	genelength=genelen;  
	ngenes=ngene;
	nmutations=0;
	nrecombinations=0;
	emergeneration=1;
	snpn=0;

	aa=0;  
	cc=0;		
	gg=0;	
	tt=0;
	GC=0;

	MLST=1;

	genes.reserve(ngenes);  
	mlstdo= pa.mlstdo;
	recsitesdo= pa.recsites;
	
	for (int i=0; i<ngenes; i++){  
		shared_ptr<Gene> geneforfor = make_shared<Gene>(genelen, pa, i);
		genes.push_back(geneforfor);
	}

	
}

void Bac::printGenes() {  
	
	for (int i=0; i<(int)genes.size(); i++) {   
		cout  << "Gene " << i << "|  "<< "Mutations: "<< genes[i]->nmutations <<  ",  Recombinations: "<< genes[i]->nrecombinations  << ",  SNPs: "<< genes[i]->snps.size() <<"  |\n" << genes[i]->getGeneText(1)  << endl;
	}
	
	
}


void Bac::mutate(Para &pa, int& snpd) { // Mutates bacterium
	
	nmutations++;
	int mutgene = ngenes*pa.getZeroone();
	
	shared_ptr<Gene> MutGene = make_shared<Gene>(*genes[mutgene]);
	MutGene->mutate(pa, snpd);
	genes[mutgene] = MutGene;
			
	snpn+=snpd;
	
	if (mlstdo>0) {
		MLST=Superpop::getMLST(genes);
	}
	
}

void Bac::recombine(int recgene, int rocus, string donatedseq, int snpchange, Para& pa) { // Recombines bacterium. Assumes non-identical recombination sequences.
	
	genes[recgene]->recombine(rocus, donatedseq);
	
	nrecombinations++; 
	
	snpn+=snpchange;
	
	if (mlstdo>0) {
		MLST=Superpop::getMLST(genes);
	}
	
	
}

void Bac::insertion(Para& pa) {
	
	int insgene = ngenes*pa.getZeroone();
	
	shared_ptr<Gene> insGene = make_shared<Gene>(*genes[insgene]);
	insGene->insertion(pa);
	genes[insgene] = insGene;

	if (mlstdo>0) {
		MLST=Superpop::getMLST(genes);
	}
	
}

void Bac::deletion(Para& pa) {
	
	int delgene = ngenes*pa.getZeroone(); 
	int len=0;
	while (len <1 || len > pa.indelmax*genes[delgene]->genstr.length()) {
			len=pa.getZipf(1,1.0);
		}
	
	Superpop::indels[-len]++;	
	
	int rocus=genes[delgene]->genstr.length()*pa.getZeroone();
	
	
	if (rocus + len > genelength) {
		len=genelength-rocus;
	}
	
	
	shared_ptr<Gene> delGene = make_shared<Gene>(*genes[delgene]);
	delGene->deletion(rocus, len);
	genes[delgene] = delGene;
	
	
	if (mlstdo>0) {
		MLST=Superpop::getMLST(genes);
	}
						
}

vector<double> Bac::getbases(int curgeneration, Para& pa) {  // Currently not in use
	
	vector<double> bases;
	bases.resize(5,0.0);
	for (int i=0; i<(int)genes.size(); i++) { 
		if (curgeneration > genes[i]->updategeneration) {
			genes[i]->updateBasecounts();
			genes[i]->updategeneration=curgeneration;
		}
		
		bases[0]+=genes[i]->aa;
		bases[1]+=genes[i]->tt;
		bases[2]+=genes[i]->gg;
		bases[3]+=genes[i]->cc;
		bases[4]+=(double)(genes[i]->gg+genes[i]->cc)/(genes[i]->aa+genes[i]->tt+genes[i]->gg+genes[i]->cc);
	}	
	bases[4]=bases[4]/(double)genes.size(); 
	
	return bases; 
}

