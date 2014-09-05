/*
 * GuideTree.h
 *
 *  Created on: Aug 28, 2014
 *      Author: root
 */

#ifndef GUIDETREE_H_
#define GUIDETREE_H_


#include "core/Definitions.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include "core/HmmException.hpp"
#include "core/BioNJ.hpp"
#include <unordered_map>
#include <vector>
#include <string>

using namespace std;

namespace EBC
{

class GuideTree
{
protected:
	Dictionary* dict;
	Sequences* inputSequences;
	unsigned int kmerSize;
	unsigned int sequenceCount;
	vector<unordered_map<string,short>*>* kmers;
	vector<double> distances;

	vector<array<unsigned int, 3> > sampledTriplets;


public:
	GuideTree(Sequences*);

	~GuideTree();

	void constructTree();

	const vector<array<unsigned int, 3> >& getSampledTriplets() const
	{
		return sampledTriplets;
	}

private:

	void extractKmers(string& seq, unordered_map<string,short>* umap);

	unsigned int commonKmerCount(unsigned int i, unsigned int j);

};

} /* namespace EBC */

#endif /* GUIDETREE_H_ */
