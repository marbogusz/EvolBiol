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
				divergenceT(0.5/*initialDivergence*/), modelParams(nullptr,nullptr, 2, 1, false, false, false, true, maths),
				fwdHmm(seq1,seq2, substModel, indelModel, Definitions::DpMatrixType::Full, nullptr,true),
				vitHmm(seq1,seq2, substModel, indelModel, Definitions::DpMatrixType::Full, nullptr,true)
{
	fwdHmm.setDivergenceTimeAndCalculateModels(divergenceT);
	vitHmm.setDivergenceTimeAndCalculateModels(divergenceT);


	forwardLnL = fwdHmm.runAlgorithm();
	vitHmm.runAlgorithm();

	modelParams.setUserDivergenceParams({divergenceT});

	//substModel->summarize();

	bfgs = new Optimizer(&modelParams,this,Definitions::OptimizationType::BFGS);

	sampleCount = Definitions::samplingPathCount;
	sampleMinCount = Definitions::samplingPathMinCount;
	sampleMaxCount = Definitions::samplingPathMaxCount;

	//FIXME - perhaps make delta dependant on the fwd lnL ???
	//e.g. forward lnl percentage ?
	sampleDeltaLnL = forwardLnL * 0.0375;//Definitions::samplingPathLnLDelta;

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

	//fwdHmm.summarize();

	//initialize
	while(sampledCount < sampleCount){
		HMMPathSample smpl;
		fwdHmm.sampleAlignment(smpl);
		tmpLnL = fwdHmm.calculateSampleLnL(smpl);

		totalSampleLnL = maths->logSum(totalSampleLnL, tmpLnL);
		//DUMP("sample : " << tmpLnL);
		if (bestLnL < tmpLnL)
			bestLnL = tmpLnL;
		if (worstLnL > tmpLnL)
			worstLnL = tmpLnL;

		samples.push_front(make_pair(tmpLnL, smpl));
		sampledCount++;
	}
	//DUMP("Best and worst");
	//DUMP(bestLnL << "\t" << worstLnL);
	//sort
	samples.sort();
	samples.unique(is_near());

	//now the head points to the lowest lnl element
	DUMP("!!!!!!SampleInitialSet sample count R1: " << sampledCount << "\tBest lnl " << bestLnL << "\tworst lnl " << worstLnL << "\ttotal lnl " << totalSampleLnL << "\tdesired lnl delta " << sampleDeltaLnL);

	//burn-in
	while(sampledCount < sampleMinCount)
	{
		HMMPathSample smpl;
		fwdHmm.sampleAlignment(smpl);
		tmpLnL = fwdHmm.calculateSampleLnL(smpl);

		if(tmpLnL > worstLnL){
			tmpPair.first = tmpLnL;
			up = std::upper_bound(samples.begin(),samples.end(), tmpPair, [](const pair<double, HMMPathSample> &val,
					const pair<double, HMMPathSample> &itp)
					{
						return ((itp.first - val.first) > 0.000001);
					});
			//move one step back
			up--;
			if (fabs(up->first - tmpLnL) > 0.000001){
				totalSampleLnL = maths->logSum(totalSampleLnL, tmpLnL);
				double tm = up->first;
				up++;
				double tu = up->first;
				samples.insert(up,make_pair(tmpLnL,smpl));
				samples.pop_front();
				worstLnL = samples.front().first;
				//cerr << "worst changes to " << worstLnL << "\t" << tm << "\t is left of upper_bound, " << tu << "\tis the bound, will insert " << tmpLnL << " left of the bound " << endl;
				//insert the element properly!

			}
			//otherwise skip!!
		}
		if(tmpLnL > bestLnL){
			bestLnL = tmpLnL;
			//cerr << "best changes to " << bestLnL << endl;
		}
		sampledCount++;
	}
	lnlDelta = bestLnL - worstLnL;

	DUMP("!!!!!!SampleInitialSet sample count R2: " << sampledCount << "\tBest lnl " << bestLnL << "\tworst lnl " << worstLnL << "\ttotal lnl " << totalSampleLnL << "\tdelta " << lnlDelta );
	//proper sampling
	while(sampledCount < sampleMaxCount && lnlDelta > sampleDeltaLnL)
	{
		//keep on sampling
		HMMPathSample smpl;
		fwdHmm.sampleAlignment(smpl);
		tmpLnL = fwdHmm.calculateSampleLnL(smpl);



		if(tmpLnL > worstLnL){
			tmpPair.first = tmpLnL;
			up = std::upper_bound(samples.begin(),samples.end(), tmpPair, [](const pair<double, HMMPathSample> &val,
					const pair<double, HMMPathSample> &itp)
					{
						return ((itp.first - val.first) > 0.000001);
					});
			//move one step back
			up--;
			if (fabs(up->first - tmpLnL) > 0.000001){
				totalSampleLnL = maths->logSum(totalSampleLnL, tmpLnL);
				double tm = up->first;
				samples.insert(++up,make_pair(tmpLnL,smpl));
				samples.pop_front();
				worstLnL = samples.front().first;
				//cerr << "worst changes to " << worstLnL << endl;
				//cerr << "worst changes to " << worstLnL << "\t" << tm << "\t is different to " << tmpLnL <<  endl;
				//insert the element properly!
			}
			//otherwise skip!!
		}
		if(tmpLnL > bestLnL){
			bestLnL = tmpLnL;
			//cerr << "best changes to " << bestLnL << endl;
		}
		sampledCount++;
		lnlDelta = bestLnL - worstLnL;
	}
	DUMP("!!!!!!SampleInitialSet sample count R3: " << sampledCount << "\tBest lnl " << bestLnL << "\tworst lnl " << worstLnL << "\ttotal lnl " << totalSampleLnL << "\tdelta " << lnlDelta );

	DUMP ("SCORED samples");
	for (auto& ent : samples){
		DUMP(ent.first << "\twith weight " << exp(ent.first - totalSampleLnL));
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
	double time = modelParams.getDivergenceTime(0);
	fwdHmm.setDivergenceTimeAndCalculateModels(time);
	/*
	for (auto &entry : samples){
		entry.first = fwdHmm.calculateSampleLnL(entry.second);
		totalSampleLnL = maths->logSum(totalSampleLnL, entry.first);
	}
	for (auto &entry : samples){
		result += (entry.first * exp(entry.first - totalSampleLnL));
	}
	 */

	//cerr << "Time " << time << "\tlnL " << result << endl;

	//return fwdHmm.calculateSampleLnL(samples.back().second) * -1.0;
	return fwdHmm.runAlgorithm();
	//return result * -1.0;
}

double PairHMMSampler::optimiseDivergenceTime()
{
	double lnl;
	//lnl = bfgs->optimize();
	//this->divergenceT = modelParams.getDivergenceTime(0);
	//return lnl;
	return runIteration();
}

} /* namespace EBC */
