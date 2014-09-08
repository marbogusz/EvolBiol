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

GuideTree::GuideTree(Sequences* is) : inputSequences(is), distMat(inputSequences->getSequenceCount())
{
	// TODO Auto-generated constructor stub
	this->dict = inputSequences->getDictionary();
	this->kmerSize = 6;
	this->sequenceCount = inputSequences->getSequenceCount();
	this->kmers = new vector<unordered_map<string,short>*>(sequenceCount);
	DEBUG("Creating the guide tree");
	this->constructTree();
}

GuideTree::~GuideTree()
{

}

void GuideTree::constructTree()
{
	unsigned int i,j;
	string currSeq;
	float identity;

	DEBUG("Extracting kmers");
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
			identity = 1.0 - commonKmerCount(i,j)/(float)(min(s1.size(),s2.size()));
			distMat.addDistance(i,j,identity);
		}

	DEBUG("initialized the DM");
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
	return commonCount;
}

} /* namespace EBC */
