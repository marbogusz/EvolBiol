/*
 * Band.cpp
 *
 *  Created on: Oct 29, 2014
 *      Author: root
 */

#include <heuristics/Band.hpp>

namespace EBC
{

Band::Band(unsigned int size) : matchBand(size), insertBand(size), deleteBand(size)
{

	posteriorLikelihoodLimit = -5;
	posteriorLikelihoodDelta = -2;
	// TODO Auto-generated constructor stub

}

Band::~Band()
{
	// TODO Auto-generated destructor stub
}

void Band::processPosteriorProbabilities(EvolutionaryPairHMM* hmm)
{
	M = hmm->getM();
	X = hmm->getX();
	Y = hmm->getY();

	//FIXME - use a defined number, same with -1000000 mins!!
	double minX;
	double minY;
	double minM;

	unsigned int xHi, xLo, yHi, yLo, mHi,mLo;

	unsigned int tmpRow;

	double tmpX, tmpY, tmpM;


	for(unsigned int col = 1; col < M->getCols(); col++)
	{
		minX = -200000;
		minY = -200000;
		minM = -200000;
		xHi=mHi=yHi=xLo=mLo=yLo=0;
		//scan the column for the smallest lnl
		for(tmpRow = 1; tmpRow < M->getRows(); tmpRow ++)
		{
			tmpX = X->getValueAt(tmpRow, col);
			if(tmpX > minX && tmpX > posteriorLikelihoodLimit)
			{
				minX = tmpX;
				xHi = tmpRow;
			}
			tmpY = Y->getValueAt(tmpRow, col);
			if(tmpY > minY && tmpY > posteriorLikelihoodLimit)
			{
				minY = tmpY;
				yHi = tmpRow;
			}
			tmpM = M->getValueAt(tmpRow, col);
			if(tmpM > minM && tmpM > posteriorLikelihoodLimit)
			{
				minM = tmpM;
				mHi = tmpRow;
			}
		}
		//got the minimums, get the range
		unsigned int tmpIdx;
		if(xHi != 0)
		{
			tmpIdx = xLo = xHi;
			while((X->getValueAt(++tmpIdx, col) - minX) > posteriorLikelihoodDelta)
				xHi = tmpIdx;
			tmpIdx = xLo;
			while((X->getValueAt(--tmpIdx, col) - minX) > posteriorLikelihoodDelta)
				xLo = tmpIdx;

		}
		if(yHi != 0)
		{
			tmpIdx = yLo = yHi;
			while((Y->getValueAt(++tmpIdx, col) - minY) > posteriorLikelihoodDelta)
				yHi = tmpIdx;
			tmpIdx = xLo;
			while((Y->getValueAt(--tmpIdx, col) - minY) > posteriorLikelihoodDelta)
				yLo = tmpIdx;

		}
		if(mHi != 0)
		{
			tmpIdx = mLo = mHi;
			while((M->getValueAt(++tmpIdx, col) - minM) > posteriorLikelihoodDelta)
				mHi = tmpIdx;
			tmpIdx = mLo;
			while((M->getValueAt(--tmpIdx, col) - minM) > posteriorLikelihoodDelta)
				mLo = tmpIdx;

		}
		this->setInsertRangeAt(col, xLo,xHi);
		this->setDeleteRangeAt(col, yLo,yHi);
		this->setMatchRangeAt(col, mLo,mHi);
	}
}

} /* namespace EBC */
