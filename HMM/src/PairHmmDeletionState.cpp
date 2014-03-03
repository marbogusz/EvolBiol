/*
 * PairHmmDeletionState.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "PairHmmDeletionState.hpp"

namespace EBC
{

PairHmmDeletionState::PairHmmDeletionState(unsigned int x, unsigned int y, double gO, double gE)
		: PairHmmState(x,y), gapOpen(gO), gapExtension(gE)
{
	initializeData();
}

PairHmmDeletionState::~PairHmmDeletionState()
{
	// TODO Auto-generated destructor stub
}

void PairHmmDeletionState::initializeData()
{
	setWholeRow(0,this->minVal);
	setWholeCol(0,this->minVal);
}


} /* namespace EBC */
