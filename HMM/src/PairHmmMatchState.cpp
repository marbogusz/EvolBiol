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
	setWholeRow(0,this->minVal);
	setWholeCol(0,this->minVal);
	//FIXME - set to some other value!
	this->matrixData[1][1].score = 0;
	//DEBUG("init state");
}


} /* namespace EBC */
