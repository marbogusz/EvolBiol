/*
 * GotohScoringMatrix.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#include "heuristics/GotohScoringMatrix.hpp"

namespace EBC
{

GotohScoringMatrix::GotohScoringMatrix(unsigned int size, double distance, Dictionary* dict) : ScoringMatrix(size, distance, dict)
{


	if(matrixSize == Definitions::nucleotideCount)
	{
		this->gapOpening = this->gapPenalty;
		this->gapExtension = gapOpening * 0.2;
	}
	else if (matrixSize == Definitions::aminoacidCount)
	{
		this->gapOpening = Definitions::blosum62gapOpening;
		this->gapExtension = Definitions::blosum62gapExtension;
	}

	DEBUG("Scoring matrix gap opening and extension scores : " << gapOpening << ", " << gapExtension << " for matrix size " << matrixSize);
}

} /* namespace EBC */
