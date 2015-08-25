/*
 * BandCalculator.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: root
 */

#include <heuristics/BandCalculator.hpp>

namespace EBC
{

BandCalculator::BandCalculator(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, SubstitutionModelBase* sm, IndelModel* im, double divergenceTime) :
		fwd(3,nullptr), seq1(s1), seq2(s2), substModel(sm), indelModel(im), time(divergenceTime)
{
	DEBUG("Band estimator running...");




	//FIXME magic numbers
	posteriorLikelihoodLimit = Definitions::bandPosteriorLikelihoodLimit;
	posteriorLikelihoodDelta = Definitions::bandPosteriorLikelihoodDelta;

	this->ptMatrix =  new PMatrixDouble(substModel);
	this->trProbs = new TransitionProbabilities(indelModel);

	bool highDivergence = false;

	//standard multipliers
	array<double,3> normalMultipliers = {{0.5,1,1.5}};
	//for high divergences
	array<double,3> highMultipliers = {{1,2,4}};


	array<double,3> &multipliers = normalMultipliers;

	band = new Band(s1->size(),s2->size());

	//default accura
	accuracy = Definitions::normalDivergenceAccuracyDelta;

	/*
	if (time > Definitions::kmerHighDivergence){
		highDivergence = true;
		multipliers = highMultipliers;
	}

	*/

	unsigned int best = 0;
	double tmpRes = std::numeric_limits<double>::max();
	double lnl;

	DUMP("Trying several forward calculations to assess the band...");
	for(unsigned int i = 0; i < fwd.size(); i++)
	{

		fwd[i] = new ForwardPairHMM(seq1,seq2, substModel,indelModel, Definitions::DpMatrixType::Full,band);
		fwd[i]->setDivergenceTimeAndCalculateModels(time*multipliers[i]);
		lnl = fwd[i]->runAlgorithm();
		DUMP("Calculation "<< i << " with divergence time " << time*multipliers[i] << " and lnL " << lnl);
		if(lnl < tmpRes)
		{
			best = i;
			tmpRes = lnl;
		}
		//bwd[i] = new BackwardPairHMM(seq1,seq2, false, substModel,indelModel, 0, Definitions::DpMatrixType::Full);
		//bwd[i]->setDivergenceTimeAndCalculateModels(time*multipliers[i]);
		//bwd[i]->runAlgorithm();
	}

	//TODO - perhaps band it as well ???
	bwd =  new BackwardPairHMM(seq1,seq2, substModel,indelModel, Definitions::DpMatrixType::Full,band);
	bwd->setDivergenceTimeAndCalculateModels(time*multipliers[best]);
	DUMP("Backward calculation runs...");
	bwd->runAlgorithm();

	//combine fwd and bwd metrics into one!

	bwd->calculatePosteriors(fwd[best]);
	this->processPosteriorProbabilities(bwd, band);

	bestTime = time*multipliers[best];
	//for high divergences likelihood surfaces tend to be flatter
	if(highDivergence){
		if(best==1){
			//divergence of 2-3
			accuracy = Definitions::highDivergenceAccuracyDelta;
		}
		else if(best==2){
			//divergence of 7-9
			accuracy = Definitions::ultraDivergenceAccuracyDelta;
		}
	}

	INFO("k-mer time: " << divergenceTime);
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

	double ccM, ccI, ccD;
	ccM = ccD = ccI = 0;
	double rowcols = M->getRows() * M->getCols();

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

		ccI += xHi-xLo+1;
		ccD += yHi-yLo+1;
		ccM += mHi-mLo+1;

		band->setDeleteRangeAt(col, yLo,yHi);
		band->setMatchRangeAt(col, mLo,mHi);

		DUMP("Match/Ins/Del bands for column " << col << "\t" << band->getMatchRangeAt(col).first <<"\t" << band->getMatchRangeAt(col).second
				<< "\t" << band->getInsertRangeAt(col).first <<"\t" << band->getInsertRangeAt(col).second
				<< "\t" << band->getDeleteRangeAt(col).first <<"\t" << band->getDeleteRangeAt(col).second);
	}
	INFO("Band promiles M:" << (int)(ccM*1000.0/rowcols) << " I:" << (int)(ccI*1000.0/rowcols) << " D:" << (int)(ccD*1000.0/rowcols));

}

double BandCalculator::getClosestDistance() {
	return this->bestTime;
}

double BandCalculator::getBrentAccuracy() {
	return this->accuracy;
}

} /* namespace EBC */


