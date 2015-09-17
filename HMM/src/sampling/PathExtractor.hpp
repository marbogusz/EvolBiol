//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================

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
