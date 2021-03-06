//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================

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
				fwdHmm(seq1,seq2, substModel, indelModel, Definitions::DpMatrixType::Full, nullptr,true),
				//vitHmm(seq1,seq2, substModel, indelModel, Definitions::DpMatrixType::Full, nullptr,true)//,
				bacHmm(seq1,seq2, substModel, indelModel, Definitions::DpMatrixType::Full, nullptr)
{
	fwdHmm.setDivergenceTimeAndCalculateModels(divergenceT);
	bacHmm.setDivergenceTimeAndCalculateModels(divergenceT);
	//vitHmm.setDivergenceTimeAndCalculateModels(divergenceT);

	fwdHmm.runAlgorithm();
	bacHmm.runAlgorithm();

	bacHmm.calculatePosteriors(&fwdHmm);

	//now we have posteriors!!!
	//Calculate MP path - requires a new type of HMM - MP hmm, right ?

	//forwardLnL = fwdHmm.runAlgorithm();
	//vitHmm.runAlgorithm();

	modelParams.setUserDivergenceParams({divergenceT});

	//substModel->summarize();

	bfgs = new Optimizer(&modelParams,this,Definitions::OptimizationType::BFGS);

	sampleCount = 10000;//Definitions::samplingPathCount;
	sampleMinCount = Definitions::samplingPathMinCount;
	sampleMaxCount = Definitions::samplingPathMaxCount;

	//FIXME - perhaps make delta dependant on the fwd lnL ???
	//e.g. forward lnl percentage ?
	sampleDeltaLnL = forwardLnL * 0.0375;//Definitions::samplingPathLnLDelta;

	totalSampleLnL = Definitions::minMatrixLikelihood;

	//sampleInitialSet();
	//remove

	//output the alignment!!!
/*
	auto val = vitHmm.getStringAlignment();

	INFO("VITERBI ALIGNMENT INITIAL");
	INFO(val.first);
	INFO(val.second);

	exit(0);
	*/
	//bacHmm.calculateMaximumPosteriorMatrix();
	//auto vb = bacHmm.getMPAlignment();

	//INFO("MP ALIGNMENT INITIAL");
	//INFO(vb.first);
	//INFO(vb.second);


	//exit(0);

	//vitHmm.getAlignment(vitSmpl);
	//bacHmm.getAlignment(vitSmpl);
	//vitHmm.getSample(seq1, seq2, vitSmpl);
	//cerr << vitHmm.calculateSampleLnL(vitSmpl) << endl;

}

PairHMMSampler::~PairHMMSampler()
{
	delete maths;
}

void PairHMMSampler::sampleInitialSet()
{
	/*
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
*/
/*
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


	//DUMP ("SCORED samples");
	//	for (auto& ent : samples){
	//		DUMP(ent.first << "\twith weight " << exp(ent.first - totalSampleLnL));
	//	}
*/
}

void PairHMMSampler::reSample()
{
	//fwdHmm.setDivergenceTimeAndCalculateModels(divergenceT);
	//fwdHmm.runAlgorithm();

	//FIXME - delete the 2 below
	//vitHmm.setDivergenceTimeAndCalculateModels(divergenceT);
	//vitHmm.runAlgorithm();

	samples.clear();
	this->sampleInitialSet();
	//vitSmpl
	//vitHmm.getAlignment(vitSmpl);
}

void PairHMMSampler::doExtraStuff(){
	/*
	fwdHmm.setDivergenceTimeAndCalculateModels(divergenceT);
	bacHmm.setDivergenceTimeAndCalculateModels(divergenceT);
	//vitHmm.setDivergenceTimeAndCalculateModels(divergenceT);

	fwdHmm.runAlgorithm();
	bacHmm.runAlgorithm();

	bacHmm.calculatePosteriors(&fwdHmm);
	bacHmm.calculateMaximumPosteriorMatrix();

	auto val = bacHmm.getMPAlignment();

	INFO("MP ALIGNMENT Final");
	INFO(val.first);
	INFO(val.second);
	*/
}

double PairHMMSampler::runIteration()
{

	double sumLnl = totalSampleLnL;
	totalSampleLnL = Definitions::minMatrixLikelihood;
	double result = Definitions::minMatrixLikelihood;;

	return 1.0;

	//double time = modelParams.getDivergenceTime(0);
	//fwdHmm.setDivergenceTimeAndCalculateModels(time);
	//fwdHmm.setDivergenceTimeAndCalculateModels(divergenceT);

	//vitHmm.setDivergenceTimeAndCalculateModels(time);
	//return fwdHmm.runAlgorithm();
/*
	for (auto &entry : samples){
		entry.first = fwdHmm.calculateSampleLnL(entry.second);
		totalSampleLnL = maths->logSum(totalSampleLnL, entry.first);
	}

	result = totalSampleLnL;
*/
/*
	for (auto &entry : samples){
		//result = maths->logSum(result, (entry.first + entry.first - totalSampleLnL));
		result = maths->logSum(result, (entry.first));
		//result = maths->logSum(result, (fwdHmm.calculateSampleLnL(entry.second)));
	}
*/
	//cerr << "Time " << time << "\tlnL " << result << endl;

	//return fwdHmm.calculateSampleLnL(samples.back().second) * -1.0;
	//vitHmm.setDivergenceTimeAndCalculateModels(time);
	//return (vitHmm.calculateSampleLnL(vitSmpl) * -1.0);
	//return fwdHmm.runAlgorithm();
	//return result * -1.0;
}

double PairHMMSampler::optimiseDivergenceTime()
{
	double lnl;
	//lnl = bfgs->optimize();
	//this->divergenceT = modelParams.getDivergenceTime(0);
	//cout << "AASDFASDGASDGA\n";
	//cout << divergenceT;
	//return lnl;

	lnl =  runIteration();

	return lnl;
}

} /* namespace EBC */
