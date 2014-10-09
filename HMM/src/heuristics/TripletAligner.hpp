/*
 * TripletAligner.hpp
 *
 *  Created on: Sep 8, 2014
 *      Author: root
 */

#ifndef TRIPLETALIGNER_HPP_
#define TRIPLETALIGNER_HPP_

#include "core/Sequences.hpp"
#include "core/DistanceMatrix.hpp"
#include "heuristics/GotohAlgorithm.hpp"
#include "core/SequenceElement.hpp"
#include <array>
#include <vector>

using namespace std;

namespace EBC
{

struct TripletAlignment
{
public:
	string seq1;
	string seq2;
	string seq3;
};

class TripletAligner
{
private:

	Sequences* inputSeqs;
	DistanceMatrix* distMat;

	unsigned int s1, s2,s3;

	array<vector<SequenceElement>, 3> triAlignment;

	void assembleFromPairs(pair<string, string>& p1,pair<string, string>& p2);

public:
	TripletAligner(Sequences* inputSeq, DistanceMatrix* dm);

	array<vector<SequenceElement>, 3> align(array<unsigned int, 3> triplet);

	array<vector<SequenceElement>, 3> align(pair<string, string>& p1,pair<string, string>& p2);
};

} /* namespace EBC */

#endif /* TRIPLETALIGNER_HPP_ */
