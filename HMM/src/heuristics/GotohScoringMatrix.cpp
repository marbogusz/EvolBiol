/*
 * GotohScoringMatrix.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#include "heuristics/GotohScoringMatrix.hpp"

namespace EBC
{

GotohScoringMatrix::GotohScoringMatrix(double distance, double gO, double gE) :
		ScoringMatrix(distance), gapOpening(gO), gapExtension(gE)
{
	// TODO Auto-generated constructor stub
}

GotohScoringMatrix::GotohScoringMatrix(double distance) : ScoringMatrix(distance)
{
	this->gapOpening = this->gapPenalty;

	//FIXME - magic number, set the gap opening to sth reasonable;
	this->gapExtension = gapOpening * 0.3;

	DEBUG("Scoring matrix gap opening and extension scores : " << gapOpening << ", " << gapExtension);
}

} /* namespace EBC */
