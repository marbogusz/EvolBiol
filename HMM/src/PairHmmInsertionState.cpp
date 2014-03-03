/*
 * PairHmmInsertionState.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "PairHmmInsertionState.hpp"

namespace EBC
{

PairHmmInsertionState::PairHmmInsertionState(unsigned int x, unsigned int y, double gO, double gE)
		: PairHmmState(x,y), gapOpen(gO), gapExtension(gE)
{
	initializeData();
}

PairHmmInsertionState::~PairHmmInsertionState()
{
	// TODO Auto-generated destructor stub
}

void PairHmmInsertionState::initializeData()
{
	setWholeRow(0,this->minVal);
	setWholeCol(0,this->minVal);
}


} /* namespace EBC */
