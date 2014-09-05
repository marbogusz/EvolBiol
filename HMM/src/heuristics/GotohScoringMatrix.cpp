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

} /* namespace EBC */
