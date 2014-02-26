/*
 * PairHmmMatchState.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "PairHmmMatchState.hpp"

namespace EBC
{

PairHmmMatchState::PairHmmMatchState(unsigned int x, unsigned int y, double gO, double gE)
		: PairHmmState(x,y), gapOpen(gO), gapExtension(gE)
{
	initializeData();
}

PairHmmMatchState::~PairHmmMatchState()
{
	// TODO Auto-generated destructor stub
}

void PairHmmMatchState::initializeData()
{
	this->matrixData[0][0].score = 1.0;
}


} /* namespace EBC */
