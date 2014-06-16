/*
 * PairwiseHmmDeletionState.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "PairwiseHmmDeleteState.hpp"
#include "DpMatrixFull.hpp";

namespace EBC
{

PairwiseHmmDeleteState::PairwiseHmmDeleteState(unsigned int x, unsigned int y)
{
    this->dpMatrix = new DpMatrixFull(x,y);
	initializeData();
}

PairwiseHmmDeleteState::PairwiseHmmDeleteState(DpMatrixBase *matrix) : dpMatrix(matrix)
{
	initializeData();
}

void PairwiseHmmDeleteState::initializeData()
{


	dpMatrix->setWholeCol(0,-100000000);
}

void PairwiseHmmDeleteState::setDirection(unsigned int i, unsigned int j)
{
	dpMatrix->setHorizontalAt(i,j);
}


} /* namespace EBC */
