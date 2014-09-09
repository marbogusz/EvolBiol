/*
 * TripletAligner.hpp
 *
 *  Created on: Sep 8, 2014
 *      Author: root
 */

#ifndef TRIPLETALIGNER_HPP_
#define TRIPLETALIGNER_HPP_

#include "core/Sequences.hpp"
#include "heuristics/GotohAlgorithm.hpp"
#include "core/SequenceElement.hpp"
#include <array>

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
	DistanceMatrix& distMat;
	unsigned int s1, s2,s3;

public:
	TripletAligner(Sequences* inputSeq, array<unsigned int, 3> triplet, DistanceMatrix& dm);

	array<vector<SequenceElement>, 3> align();
};

} /* namespace EBC */

#endif /* TRIPLETALIGNER_HPP_ */
