/*
 * PairwiseHmmInsertionState.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "PairwiseHmmInsertState.hpp"
#include "DpMatrixFull.hpp"

namespace EBC
{

PairwiseHmmInsertState::PairwiseHmmInsertState(unsigned int x, unsigned int y)
{
    this->dpMatrix = new DpMatrixFull(x,y);
	initializeData();
}

PairwiseHmmInsertState::PairwiseHmmInsertState(DpMatrixBase *matrix)
{
	this->dpMatrix = matrix;
	initializeData();
}

void PairwiseHmmInsertState::initializeData()
{
	dpMatrix->setWholeRow(0,-100000);
}

void PairwiseHmmInsertState::setDirection(unsigned int i, unsigned int j)
{
	dpMatrix->setVerticalAt(i,j);
}

PairwiseHmmInsertState::~PairwiseHmmInsertState()
{
}

} /* namespace EBC */
