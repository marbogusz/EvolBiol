/*
 * GotohScoringMatrix.h
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#ifndef GOTOHSCORINGMATRIX_H_
#define GOTOHSCORINGMATRIX_H_

#include "heuristics/ScoringMatrix.hpp"
#include <string>

using namespace std;

namespace EBC
{

class GotohScoringMatrix: public EBC::ScoringMatrix
{
protected:

	double gapOpening;
	double gapExtension;

	//Pairwise alignment references
	string seq1;
	string seq2;

public:
	GotohScoringMatrix(double disance, double gapOpening, double gapExtension);

	double getScoreByIndex(unsigned int i, unsigned int j)
	{
		return this->getScore(seq1[i], seq2[j]);
	}

	void setStringPair(string& s1, string& s2)
	{
		seq1 = s1;
		seq2 = s2;
	}

	double getGapExtension() const
	{
		return gapExtension;
	}

	void setGapExtension(double gapExtension)
	{
		this->gapExtension = gapExtension;
	}

	double getGapOpening() const
	{
		return gapOpening;
	}

	void setGapOpening(double gapOpening)
	{
		this->gapOpening = gapOpening;
	}
};

} /* namespace EBC */
#endif /* GOTOHSCORINGMATRIX_H_ */
