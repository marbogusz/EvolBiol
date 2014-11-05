/*
 * GuideTree.cpp
 *
 *  Created on: Aug 28, 2014
 *      Author: root
 */

#include "heuristics/GuideTree.hpp"
#include "heuristics/TripletSamplingTree.hpp"

namespace EBC
{

GuideTree::GuideTree(Sequences* is) : inputSequences(is)
{
	// TODO Auto-generated constructor stub
	distMat = new DistanceMatrix(inputSequences->getSequenceCount());
	this->dict = inputSequences->getDictionary();
	if(dict->getAlphabetSize() == Definitions::nucleotideCount)
		this->kmerSize = 7;
	else
		this->kmerSize = 4;
	this->sequenceCount = inputSequences->getSequenceCount();
	this->kmers = new vector<unordered_map<string,short>*>(sequenceCount);
	DEBUG("Creating guide tree");
	this->constructTree();
}

GuideTree::~GuideTree()
{
	delete kmers;
}

double GuideTree::kimuraDist(double id)
{
		double p = 1 - id;
	// Typical case: use Kimura's empirical formula
		//if (p < 0.75)
			return -log(1 - p - (p*p)/5);

	// Per ClustalW, return 10.0 for anything over 93%
	//	if (p > 0.93)
	//		return 10.0;

	// If p >= 0.75, use table lookup
	//	assert(p <= 1 && p >= 0.75);
	// Thanks for Michael Hoel for pointing out a bug
	// in the table index calculation in versions <= 3.52.
	//	int iTableIndex = (int) ((p - 0.75)*1000 + 0.5);
	//	if (iTableIndex < 0 || iTableIndex >= iTableEntries)
	//		Quit("Internal error in MSADistKimura::ComputeDist");

	//	return dayhoff_pams[iTableIndex] / 100.0;
}

void GuideTree::constructTree()
{
	unsigned int i,j;
	string currSeq;
	double identity, estIdentity;



	DEBUG("Extracting k-mers");
	for(i = 0; i< sequenceCount; i++)
	{
		(*kmers)[i] = new unordered_map<string,short>();
		currSeq = inputSequences->getRawSequenceAt(i);
		extractKmers(currSeq, (*kmers)[i]);
	}
	for(i = 0; i< sequenceCount; i++)
		for(j = i+1; j< sequenceCount; j++)
		{
			string s1 = inputSequences->getRawSequenceAt(i);
			string s2 = inputSequences->getRawSequenceAt(j);
			identity = 1.0 - commonKmerCount(i,j)/((double)(min(s1.size(),s2.size())));
			//estIdentity = log(0.02 + identity)/4.12 + 0.995;
			//kimura = kimuraDist(estIdentity);

			//DEBUG("k-mer distance between seq. " << i << " and " << j << " is " << identity << " " << -log(0.02 + identity) << " " << -log(0.1 + identity)  << " "<< estIdentity << " " << kimura );

			//if(dict->getAlphabetSize() == Definitions::nucleotideCount)
				estIdentity = identity;
				//estIdentity = nucFunction(identity);
			//else if(dict->getAlphabetSize() == Definitions::aminoacidCount)
			//	estIdentity = aaFunction(identity);

			DEBUG("k-mer distance between seq. " << i << " and " << j << " is " << identity << " " << estIdentity );

			distMat->addDistance(i,j,estIdentity);
		}

	DEBUG("Initialized distance matrix");
	BioNJ nj(sequenceCount, distMat);
	newickTree = nj.calculate();

	/*
	TripletSamplingTree tst(distMat);
	tst.fromNewick(tree);
	DEBUG("Created the tree from Newick");

	//FIXME - do it nicely - guide tree returns a newick string, sampling tree build based on gt
	sampledTriplets = tst.sampleFromDM();
	DEBUG("Obtained triplets");
	 */
}

void GuideTree::extractKmers(string& seq, unordered_map<string, short>* umap)
{
	string kmer;



	for(unsigned int i = 0; i< seq.size() - kmerSize; i++)
	{
		kmer = seq.substr(i, kmerSize);
		++((*umap)[kmer]);
	}
}

unsigned int GuideTree::commonKmerCount(unsigned int i, unsigned int j)
{
	unordered_map<string, short>* m1 = (*kmers)[i];
	unordered_map<string, short>* m2 = (*kmers)[j];
	unsigned int commonCount = 0;

	for(auto it = m1->begin(); it != m1->end(); it++)
	{
		commonCount += std::min((short)((*m2)[it->first]), (short)(it->second));
	}
	//DEBUG("Common k-mer count between seq. " << i << " and " << j << " is " << commonCount);
	return commonCount;


}

} /* namespace EBC */
