/*
 * GuideTree.h
 *
 *  Created on: Aug 28, 2014
 *      Author: root
 */

#ifndef GUIDETREE_H_
#define GUIDETREE_H_


#include "Definitions.hpp"
#include "Maths.hpp"
#include "Dictionary.hpp"
#include "Sequences.hpp"
#include "HmmException.hpp"
#include "BioNJ.hpp"
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

public:
	GuideTree(Sequences*);

	~GuideTree();

	void constructTree();

private:

	void extractKmers(string& seq, unordered_map<string,short>* umap);

	unsigned int commonKmerCount(unsigned int i, unsigned int j);

};

} /* namespace EBC */

#endif /* GUIDETREE_H_ */
