/*
 * PairwiseHmmMatchState.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "PairwiseHmmMatchState.hpp"
#include "DpMatrixFull.hpp"

namespace EBC
{

PairwiseHmmMatchState::PairwiseHmmMatchState(unsigned int x, unsigned int y)
{
    this->dpMatrix = new DpMatrixFull(x,y);
	initializeData();
}

PairwiseHmmMatchState::PairwiseHmmMatchState(DpMatrixBase *matrix)
{
	this->dpMatrix = matrix;
	initializeData();
}

PairwiseHmmMatchState::~PairwiseHmmMatchState()
{
}

void PairwiseHmmMatchState::initializeData()
{
	dpMatrix->setWholeCol(0,-100000);
	dpMatrix->setWholeRow(0,-100000);
	dpMatrix->setValue(0,0,0);
}

void PairwiseHmmMatchState::setDirection(unsigned int i, unsigned int j)
{
	this->dpMatrix->setDiagonalAt(i,j);
}

} /* namespace EBC */
