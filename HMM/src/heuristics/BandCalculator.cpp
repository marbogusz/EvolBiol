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
		fwd(3,nullptr), seq1(s1), seq2(s2), substModel(sm), indelModel(im), time(divergenceTime)
{
	// TODO Auto-generated constructor stub
	//FIXME magic numbers
	posteriorLikelihoodLimit = -3;
	posteriorLikelihoodDelta = -9;

	this->ptMatrix =  new PMatrixDouble(substModel);
	this->trProbs = new TransitionProbabilities(indelModel);
	//FIXME - magic numbers
	array<double,3> multipliers = {{0.5,1,1.5}};
	unsigned int best = 0;
	double tmpRes = std::numeric_limits<double>::max();
	double lnl;

	for(unsigned int i = 0; i < fwd.size(); i++)
	{

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

	DEBUG("Best " << best << " time " << time*multipliers[best]);

	bwd =  new BackwardPairHMM(seq1,seq2, substModel,indelModel, Definitions::DpMatrixType::Full);
	DEBUG("BWD Set time");
	bwd->setDivergenceTime(time*multipliers[best]);
	DEBUG("BWD Run");
	bwd->runAlgorithm();

	//combine fwd and bwd metrics into one!
	band = new Band(s2.size()+1);

	FileLogger::DebugLogger() << "Posterior probabilities:\n";

	bwd->calculatePosteriors(fwd[best]);
	this->processPosteriorProbabilities(bwd, band);
}

BandCalculator::~BandCalculator()
{
	delete bwd;

	for(int i=0; i< fwd.size(); i++)
	{
		delete fwd[i];
	}

	delete trProbs;
	delete ptMatrix;

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

	//cumulative posterior likelihood
	double cpl = posteriorLikelihoodLimit + posteriorLikelihoodDelta;

	int xHi, xLo, yHi, yLo, mHi,mLo, rowCount;

	unsigned int tmpRow;

	double tmpX, tmpY, tmpM;

	rowCount = M->getRows();

	for(unsigned int col = 0; col < M->getCols(); col++)
	{
		xHi=mHi=yHi=xLo=mLo=yLo=-1;
		tmpRow = 0;
		while(tmpRow < rowCount && X->getValueAt(tmpRow, col) < cpl )
			tmpRow++;
		if (tmpRow != rowCount)
		{
			//found a value
			xLo = tmpRow;
			tmpRow = rowCount-1;
			while(tmpRow >= 0 && X->getValueAt(tmpRow, col) < cpl)
				tmpRow--;
			if (tmpRow > 0)
				xHi = tmpRow;
		}

		tmpRow = 0;
		while(tmpRow < rowCount && Y->getValueAt(tmpRow, col) < cpl)
			tmpRow++;
		if (tmpRow != rowCount)
		{
			//found a value
			yLo = tmpRow;
			tmpRow = rowCount-1;
			while(tmpRow >= 0 && Y->getValueAt(tmpRow, col) < cpl)
				tmpRow--;
			if (tmpRow > 0)
				yHi = tmpRow;
		}

		tmpRow = 0;
		while(tmpRow < rowCount && M->getValueAt(tmpRow, col) < cpl)
			tmpRow++;
		if (tmpRow != rowCount)
		{
			//found a value
			mLo = tmpRow;
			tmpRow = rowCount-1;
			while(tmpRow >= 0 && M->getValueAt(tmpRow, col) < cpl)
				tmpRow--;
			if (tmpRow > 0)
				mHi = tmpRow;
		}



		band->setInsertRangeAt(col, xLo,xHi);
		band->setDeleteRangeAt(col, yLo,yHi);
		band->setMatchRangeAt(col, mLo,mHi);

		DEBUG("M/X/Y bands " << col << "\t" << band->getMatchRangeAt(col).first <<"\t" << band->getMatchRangeAt(col).second
				<< "\t" << band->getInsertRangeAt(col).first <<"\t" << band->getInsertRangeAt(col).second
				<< "\t" << band->getDeleteRangeAt(col).first <<"\t" << band->getDeleteRangeAt(col).second);
	}
}


/*
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

	int xHi, xLo, yHi, yLo, mHi,mLo, rowCount;

	unsigned int tmpRow;

	double tmpX, tmpY, tmpM;

	rowCount = M->getRows();

	for(unsigned int col = 0; col < M->getCols(); col++)
	{
		minX = Definitions::minMatrixLikelihood *2;
		minY = Definitions::minMatrixLikelihood *2;
		minM = Definitions::minMatrixLikelihood *2;
		xHi=mHi=yHi=xLo=mLo=yLo=-1;
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
		if(xHi != -1)
		{
			tmpIdx = xLo = xHi;
			while(++tmpIdx < rowCount && (X->getValueAt(tmpIdx, col) - minX) > posteriorLikelihoodDelta)
				xHi = tmpIdx;
			tmpIdx = xLo;
			while(--tmpIdx >=0 && (X->getValueAt(tmpIdx, col) - minX) > posteriorLikelihoodDelta)
				xLo = tmpIdx;

		}
		if(yHi != -1)
		{
			tmpIdx = yLo = yHi;
			while(++tmpIdx < rowCount && (Y->getValueAt(tmpIdx, col) - minY) > posteriorLikelihoodDelta)
				yHi = tmpIdx;
			tmpIdx = xLo;
			while(--tmpIdx >=0 && (Y->getValueAt(tmpIdx, col) - minY) > posteriorLikelihoodDelta)
				yLo = tmpIdx;

		}
		if(mHi != -1)
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

		//int min = std::min(std::min(xLo < 0 ? 1 : xLo,yLo < 0 ? 1 : yLo),mLo < 0 ? 1 : mLo);
		int max = std::max(std::max(xHi,yHi),mHi);


		int Xmin = std::min(std::min(xLo < 1 ? 10000 : xLo,yLo < 10000 ? 100000 : yLo),mLo < 1 ? 10000 : mLo);


		int Ymin = std::min(std::min(xLo < 0 ? 0 : xLo,yLo < 0 ? 0 : yLo),mLo < 0 ? 0 : mLo);


		int Mmin = std::min(std::min(xLo < 1 ? 10000 : xLo,yLo < 1 ? 10000 : yLo),mLo < 1 ? 10000 : mLo);




		//FIXME - do some proper estimation instead doing band intersection
		//band->setInsertRangeAt(col, Xmin, max);
		//band->setDeleteRangeAt(col, Ymin,max);
		//band->setMatchRangeAt(col, Mmin,max);
		//band->setInsertRangeAt(col, 1,rowCount-1);
		//band->setDeleteRangeAt(col, 0,rowCount-1);
		//band->setMatchRangeAt(col, 1,rowCount-1);
		DEBUG("M/X/Y bands " << col << "\t" << band->getMatchRangeAt(col).first <<"\t" << band->getMatchRangeAt(col).second
				<< "\t" << band->getInsertRangeAt(col).first <<"\t" << band->getInsertRangeAt(col).second
				<< "\t" << band->getDeleteRangeAt(col).first <<"\t" << band->getDeleteRangeAt(col).second);
	}
}
*/

} /* namespace EBC */
