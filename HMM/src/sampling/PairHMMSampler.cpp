/*
 * PairHMMSampler.cpp
 *
 *  Created on: Jan 26, 2015
 *      Author: marcin
 */

#include "sampling/PairHMMSampler.hpp"

namespace EBC {

PairHMMSampler::PairHMMSampler(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2,
		SubstitutionModelBase* smdl, IndelModel* imdl, double initialDivergence) :
				maths(new Maths()), seq1(s1), seq2(s2), substModel(smdl), indelModel(imdl),
				divergenceT(initialDivergence), modelParams(nullptr,nullptr, 2, 1, false, false, false, true, maths),
				bfgs(&modelParams,this,Definitions::OptimizationType::BFGS),
				fwdHmm(seq1,seq2, substModel, indelModel, Definitions::DpMatrixType::Full, nullptr,true)

{
	fwdHmm.setDivergenceTimeAndCalculateModels(divergenceT);
	forwardLnL = fwdHmm.runAlgorithm();

	sampleCount = Definitions::samplingPathCount;
	sampleMinCount = Definitions::samplingPathMinCount;
	sampleMaxCount = Definitions::samplingPathMaxCount;

	//FIXME - perhaps make delta dependant on the fwd lnL ???
	//e.g. forward lnl percentage ?
	sampleDeltaLnL = Definitions::samplingPathLnLDelta;

	totalSampleLnL = Definitions::minMatrixLikelihood;

	sampleInitialSet();

}

PairHMMSampler::~PairHMMSampler()
{
	delete maths;
}

void PairHMMSampler::sampleInitialSet()
{
	double bestLnL = Definitions::minMatrixLikelihood;
	double worstLnL = 0;
	double lnlDelta;
	double tmpLnL;
	HMMPathSample dummy;
	pair<double, HMMPathSample> tmpPair = make_pair(bestLnL, dummy);

	list<pair<double, HMMPathSample> >::iterator up;

	unsigned int sampledCount = 0;

	fwdHmm.summarize();

	//initialize
	while(sampledCount < sampleCount){
		HMMPathSample smpl;
		fwdHmm.sampleAlignment(smpl);
		tmpLnL = fwdHmm.calculateSampleLnL(smpl);

		totalSampleLnL = maths->logSum(totalSampleLnL, tmpLnL);

		if (bestLnL < tmpLnL)
			bestLnL = tmpLnL;
		if (worstLnL > tmpLnL)
			worstLnL = tmpLnL;

		samples.push_front(make_pair(tmpLnL, smpl));
		sampledCount++;
	}
	//sort
	samples.sort();
	//now the head points to the lowest lnl element

	//burn-in
	while(sampledCount < sampleMinCount)
	{
		HMMPathSample smpl;
		fwdHmm.sampleAlignment(smpl);
		tmpLnL = fwdHmm.calculateSampleLnL(smpl);

		totalSampleLnL = maths->logSum(totalSampleLnL, tmpLnL);

		if(tmpLnL > worstLnL){
			samples.pop_front();
			worstLnL = tmpLnL;
			//insert the element properly!
			tmpPair.first = tmpLnL;
			up = std::upper_bound(samples.begin(),samples.end(), tmpPair);
			samples.insert(up,make_pair(tmpLnL,smpl));
		}
		if(tmpLnL > bestLnL){
			bestLnL = tmpLnL;
		}
		sampledCount++;
	}
	lnlDelta = bestLnL - worstLnL;
	//proper sampling
	while(sampledCount < sampleMinCount && lnlDelta > sampleDeltaLnL)
	{
		//keep on sampling
		HMMPathSample smpl;
		fwdHmm.sampleAlignment(smpl);
		tmpLnL = fwdHmm.calculateSampleLnL(smpl);

		totalSampleLnL = maths->logSum(totalSampleLnL, tmpLnL);

		if(tmpLnL > worstLnL){
			samples.pop_front();
			worstLnL = tmpLnL;
			//insert the element properly!

			tmpPair.first = tmpLnL;
			up = std::upper_bound(samples.begin(),samples.end(), tmpPair);
			samples.insert(up,make_pair(tmpLnL,smpl));
		}
		if(tmpLnL > bestLnL){
			bestLnL = tmpLnL;
		}
		sampledCount++;
		lnlDelta = bestLnL - worstLnL;
	}
}

void PairHMMSampler::reSample()
{
}

double PairHMMSampler::runIteration()
{
	double sumLnl = totalSampleLnL;
	totalSampleLnL = Definitions::minMatrixLikelihood;
	double result = 0;
	fwdHmm.setDivergenceTimeAndCalculateModels(modelParams.getDivergenceTime(0));
	for (auto entry : samples){
		entry.first = fwdHmm.calculateSampleLnL(entry.second);
		totalSampleLnL = maths->logSum(totalSampleLnL, entry.first);
		result += (entry.first * exp(entry.first - sumLnl));
	}
	return result;
}

double PairHMMSampler::optimiseDivergenceTime()
{
	double lnl;
	lnl = bfgs.optimize();
	this->divergenceT = modelParams.getDivergenceTime(0);
	return lnl;
}

} /* namespace EBC */
