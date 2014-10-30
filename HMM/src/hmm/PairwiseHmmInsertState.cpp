/*
 * PairwiseHmmInsertionState.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "hmm/PairwiseHmmInsertState.hpp"
#include "hmm/DpMatrixFull.hpp"

namespace EBC
{

PairwiseHmmInsertState::PairwiseHmmInsertState(unsigned int x, unsigned int y)
{
	this->rows =x;
	this->cols =y;
    this->dpMatrix = new DpMatrixFull(x,y);
	//initializeData();
}

PairwiseHmmInsertState::PairwiseHmmInsertState(DpMatrixBase *matrix)
{
	this->dpMatrix = matrix;
	//initializeData();
}

void PairwiseHmmInsertState::initializeData(bool backwards)
{
	/*
	if(!backwards)
		dpMatrix->setWholeRow(0,-100000);
	else
	{
		dpMatrix->setWholeCol(this->cols-1,-100000);
		dpMatrix->setWholeRow(this->rows-1,-100000);
		dpMatrix->setValue(rows-1,cols-1,0);
	}
	*/
}

void PairwiseHmmInsertState::setDirection(unsigned int i, unsigned int j)
{
	dpMatrix->setVerticalAt(i,j);
}

PairwiseHmmInsertState::~PairwiseHmmInsertState()
{
}

} /* namespace EBC */
