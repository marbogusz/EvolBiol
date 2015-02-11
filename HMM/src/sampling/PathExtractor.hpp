/*
 * PathExtractor.hpp
 *
 *  Created on: Feb 11, 2015
 *      Author: marcin
 */

#ifndef SAMPLING_PATHEXTRACTOR_HPP_
#define SAMPLING_PATHEXTRACTOR_HPP_

#include "core/SequenceElement.hpp"
#include "core/Maths.hpp"
#include "hmm/EvolutionaryPairHMM.hpp"

#include <string>
#include <vector>
#include <iostream>

using namespace std;

namespace EBC {

class PathExtractor {
	//reference sequences
	string s1;
	string s2;

	vector<SequenceElement*> v1;
	vector<SequenceElement*> v2;
	//lengths
	unsigned int l1;
	unsigned int l2;

	vector<pair<string,string> > sPaths;

	vector<pair<vector<SequenceElement*>,vector<SequenceElement*> > > vPaths;

public:
	PathExtractor(string seq1, string seq2, vector<SequenceElement*>* vp1, vector<SequenceElement*>* vp2);

	void genAllSeqs(string c1, string c2, unsigned int i1, unsigned int i2, unsigned int s);

	void genAllSeqs(vector<SequenceElement*> c1, vector<SequenceElement*> c2, SequenceElement* gap, unsigned int i1, unsigned int i2, unsigned int s);

	double getSumLnl(EvolutionaryPairHMM* phmm, Maths* mt);

	virtual ~PathExtractor();
};

} /* namespace EBC */

#endif /* SAMPLING_PATHEXTRACTOR_HPP_ */
