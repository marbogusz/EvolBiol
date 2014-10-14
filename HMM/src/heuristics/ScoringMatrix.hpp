/*
 * ScoringMatrix.h
 *
 *  Created on: Sep 25, 2013
 *      Author: mbogusz
 */

#ifndef SCORINGMATRIX_H_
#define SCORINGMATRIX_H_


//Needleman Wunsch scoring matrix
//TODO - legacy from a previous project - clean up by having 1 matrix class for Gotoh
#include <string>
#include "core/Definitions.hpp"
#include "core/Dictionary.hpp"

using namespace std;

namespace EBC
{

class ScoringMatrix
{
protected:

	double gapPenalty;
	double** scores;

	Dictionary* dict;

	unsigned int matrixSize;

	void scoresFromBLOSUM();

	void scoresFromDistanceJC(double distance);

public:
	ScoringMatrix();

	~ScoringMatrix();

	ScoringMatrix(unsigned int size, double distance,Dictionary * dict);

	double getScore(char& a, char& b);

	double getGapPenalty() const
	{
		return gapPenalty;
	}

	void setGapPenalty(double gapPenalty)
	{
		this->gapPenalty = gapPenalty;
	}

	double getAlignmentScore(pair<string, string>* alignment);
};

} /* namespace EBC */
#endif /* SCORINGMATRIX_H_ */
