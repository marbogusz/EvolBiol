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
	//setWholeRow(0,this->minVal);
	//setWholeCol(0,this->minVal);
	double logOpen = log(gapOpen);
	double logExtension = log(gapExtension);



	for (int i=0;i<xSize; i++)
	{
		this->matrixData[i][0].score = logOpen + i*logExtension;
	}

	for (int j=0;j<ySize; j++)
	{
		this->matrixData[0][j].score = logOpen + j*logExtension;
	}

}


} /* namespace EBC */
