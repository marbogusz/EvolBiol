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
#include "core/SequenceElement.hpp"
#include <array>
#include <vector>
#include <cmath>

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

	double posteriorTsh;
	bool usePosteriors;

	unsigned int s1, s2,s3;

	array<vector<SequenceElement>, 3> triAlignment;

	void assembleFromPairs(pair<string, string>& p1,pair<string, string>& p2);

public:
	TripletAligner(Sequences* inputSeq, DistanceMatrix* dm, double postTsh);

	//array<vector<SequenceElement>, 3> align(array<unsigned int, 3> triplet);

	//array<vector<SequenceElement>, 3> align(pair<string, string>& p1,pair<string, string>& p2);

	array<vector<SequenceElement*>*, 3>* align(pair<vector<SequenceElement*>*, vector<SequenceElement*>* >* p1, pair<vector<SequenceElement*>*, vector<SequenceElement*>* >* p2);

	array<vector<unsigned char>*, 3>* align(pair<vector<unsigned char>*, vector<unsigned char>* >* p1, pair<vector<unsigned char>*, vector<unsigned char>* >* p2);

	array<vector<unsigned char>*, 3> alignPosteriors(pair<vector<unsigned char>*, vector<unsigned char>* > p1,
			pair<vector<unsigned char>*, vector<unsigned char>* > p2, vector<double>* postP1, vector<double>* postP2);
};

} /* namespace EBC */

#endif /* TRIPLETALIGNER_HPP_ */
