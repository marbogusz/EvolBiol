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
	double logOpen = log(gapOpen);
	double logExtension = log(gapExtension);



	/*for (int i=0;i<xSize; i++)
	{
		this->matrixData[i][0].score = logOpen + i*logExtension;
	}

	for (int j=0;j<ySize; j++)
	{
		this->matrixData[0][j].score = logOpen + j*logExtension;
	}
*/
	this->setWholeCol(0,-100000000);
	this->setWholeRow(0,-100000000);
	this->matrixData[0][0].score = 0;
	//DEBUG("init state");
}

void PairHmmMatchState::setDirection(unsigned int i, unsigned int j)
{
	this->matrixData[i][j].diag = true;
}


} /* namespace EBC */
