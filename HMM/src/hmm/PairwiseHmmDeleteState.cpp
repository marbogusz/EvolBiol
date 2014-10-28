/*
 * PairwiseHmmDeletionState.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "hmm/PairwiseHmmDeleteState.hpp"
#include "hmm/DpMatrixFull.hpp"

namespace EBC
{

PairwiseHmmDeleteState::PairwiseHmmDeleteState(unsigned int x, unsigned int y)
{
	this->rows =x;
	this->cols =y;
	this->dpMatrix = new DpMatrixFull(x,y);
	//initializeData();
}

PairwiseHmmDeleteState::PairwiseHmmDeleteState(DpMatrixBase *matrix)
{
	this->dpMatrix = matrix;
	//initializeData();
}

void PairwiseHmmDeleteState::initializeData(bool backwards)
{
	if(!backwards)
		dpMatrix->setWholeCol(0,-100000);
	else
	{
		dpMatrix->setWholeCol(this->cols-1,-100000);
		dpMatrix->setWholeRow(this->rows-1,-100000);
		dpMatrix->setValue(rows-1,cols-1,0);
	}
}

void PairwiseHmmDeleteState::setDirection(unsigned int i, unsigned int j)
{
	dpMatrix->setHorizontalAt(i,j);
}

PairwiseHmmDeleteState::~PairwiseHmmDeleteState()
{
}

} /* namespace EBC */
