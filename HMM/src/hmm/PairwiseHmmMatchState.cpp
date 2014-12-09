/*
 * PairwiseHmmMatchState.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "hmm/PairwiseHmmMatchState.hpp"
#include "hmm/DpMatrixFull.hpp"

namespace EBC
{

PairwiseHmmMatchState::PairwiseHmmMatchState(unsigned int x, unsigned int y)
{
	this->rows =x;
	this->cols =y;
    this->dpMatrix = new DpMatrixFull(x,y);
    stateId = Definitions::StateId::Match;
	//initializeData();
}

PairwiseHmmMatchState::PairwiseHmmMatchState(DpMatrixBase *matrix)
{
	this->dpMatrix = matrix;
	//initializeData();
}

PairwiseHmmMatchState::~PairwiseHmmMatchState()
{
}

void PairwiseHmmMatchState::initializeData(bool backwards)
{
	if (!backwards)
	{
		//dpMatrix->setWholeCol(0,-100000);
		//dpMatrix->setWholeRow(0,-100000);
		dpMatrix->setValue(0,0,0);
	}
	/*else
	{
		dpMatrix->setWholeCol(this->cols-1,-100000);
		dpMatrix->setWholeRow(this->rows-1,-100000);
		dpMatrix->setValue(rows-1,cols-1,0);
	}
	*/
}

void PairwiseHmmMatchState::setDirection(unsigned int i, unsigned int j)
{
	this->dpMatrix->setDiagonalAt(i,j);
}

} /* namespace EBC */
