/*
 * BandCalculator.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: root
 */

#include <heuristics/BandCalculator.hpp>

namespace EBC
{

BandCalculator::BandCalculator(vector<SequenceElement>& s1, vector<SequenceElement>& s2, SubstitutionModelBase* sm, IndelModel* im, double divergenceTime) :
		fwd(3), seq1(s1), seq2(s2), substModel(sm), indelModel(im), time(divergenceTime)
{
	// TODO Auto-generated constructor stub
	//FIXME magic numbers
	posteriorLikelihoodLimit = -5;
	posteriorLikelihoodDelta = -2;

	this->ptMatrix =  new PMatrixDouble(substModel);
	this->trProbs = new TransitionProbabilities(indelModel);
	//FIXME - magic numbers
	array<double,3> multipliers = {{0.5,1,1.5}};
	unsigned int best = 0;
	double tmpRes = std::numeric_limits<double>::max();
	double lnl;

	for(unsigned int i = 0; i < fwd.size(); i++)
	{

		DEBUG(i << " with divergence time " << divergenceTime);
		fwd[i] = new ForwardPairHMM(seq1,seq2, substModel,indelModel, Definitions::DpMatrixType::Full);
		fwd[i]->setDivergenceTime(time*multipliers[i]);
		lnl = fwd[i]->runAlgorithm();
		DEBUG(i << " with divergence time " << divergenceTime << " lnl " << lnl);
		if(lnl < tmpRes)
		{
			best = i;
			tmpRes = lnl;
		}
		//bwd[i] = new BackwardPairHMM(seq1,seq2, false, substModel,indelModel, 0, Definitions::DpMatrixType::Full);
		//bwd[i]->setDivergenceTime(time*multipliers[i]);
		//bwd[i]->runAlgorithm();
	}

	DEBUG("Best " << best);

	bwd =  new BackwardPairHMM(seq1,seq2, substModel,indelModel, Definitions::DpMatrixType::Full);
	DEBUG("BWD Set time");
	bwd->setDivergenceTime(time*multipliers[best]);
	DEBUG("BWD Run");
	bwd->runAlgorithm();

	//combine fwd and bwd metrics into one!
	band = new Band(s2.size());
	bwd->calculatePosteriors(fwd[best]);
	this->processPosteriorProbabilities(bwd, band);
}

BandCalculator::~BandCalculator()
{
	for(int i=0; i< fwd.size(); i++)
	{
		delete fwd[i];
	}
	delete bwd;
}

void BandCalculator::processPosteriorProbabilities(BackwardPairHMM* hmm, Band* band)
{
	//Match state
	PairwiseHmmStateBase* M;
	//Insert state
	PairwiseHmmStateBase* X;
	//Delete state
	PairwiseHmmStateBase* Y;

	M = hmm->getM();
	X = hmm->getX();
	Y = hmm->getY();

	//FIXME - use a defined number, same with -1000000 mins!!
	double minX;
	double minY;
	double minM;

	unsigned int xHi, xLo, yHi, yLo, mHi,mLo, rowCount;

	unsigned int tmpRow;

	double tmpX, tmpY, tmpM;

	rowCount = M->getRows();

	for(unsigned int col = 0; col < M->getCols(); col++)
	{
		minX = -200000;
		minY = -200000;
		minM = -200000;
		xHi=mHi=yHi=xLo=mLo=yLo=0;
		//scan the column for the smallest lnl
		for(tmpRow = 0; tmpRow < rowCount; tmpRow ++)
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
		int tmpIdx;
		if(xHi != 0)
		{
			tmpIdx = xLo = xHi;
			while(++tmpIdx < rowCount && (X->getValueAt(tmpIdx, col) - minX) > posteriorLikelihoodDelta)
				xHi = tmpIdx;
			tmpIdx = xLo;
			while(--tmpIdx >=0 && (X->getValueAt(tmpIdx, col) - minX) > posteriorLikelihoodDelta)
				xLo = tmpIdx;

		}
		if(yHi != 0)
		{
			tmpIdx = yLo = yHi;
			while(++tmpIdx < rowCount && (Y->getValueAt(tmpIdx, col) - minY) > posteriorLikelihoodDelta)
				yHi = tmpIdx;
			tmpIdx = xLo;
			while(--tmpIdx >=0 && (Y->getValueAt(tmpIdx, col) - minY) > posteriorLikelihoodDelta)
				yLo = tmpIdx;

		}
		if(mHi != 0)
		{
			tmpIdx = mLo = mHi;
			while(++tmpIdx < rowCount && (M->getValueAt(tmpIdx, col) - minM) > posteriorLikelihoodDelta)
				mHi = tmpIdx;
			tmpIdx = mLo;
			while(--tmpIdx >=0 && (M->getValueAt(tmpIdx, col) - minM) > posteriorLikelihoodDelta)
				mLo = tmpIdx;

		}
		band->setInsertRangeAt(col, xLo,xHi);
		band->setDeleteRangeAt(col, yLo,yHi);
		band->setMatchRangeAt(col, mLo,mHi);
		DEBUG("M band " << col << "\t" << mLo <<"\t" << mHi);
	}
}


} /* namespace EBC */
