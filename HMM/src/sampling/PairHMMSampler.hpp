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
 * PairHMMSampler.hpp
 *
 *  Created on: Jan 26, 2015
 *      Author: marcin
 */

#ifndef SAMPLING_PAIRHMMSAMPLER_HPP_
#define SAMPLING_PAIRHMMSAMPLER_HPP_

#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include "core/HmmException.hpp"
#include "core/IOptimizable.hpp"
#include "core/Optimizer.hpp"

#include "models/SubstitutionModelBase.hpp"
#include "models/NegativeBinomialGapModel.hpp"

#include "hmm/ForwardPairHMM.hpp"
#include "hmm/ViterbiPairHMM.hpp"
#include "hmm/BackwardPairHMM.hpp"

#include "sampling/HMMPathSample.hpp"

#include <list>
#include <functional>

namespace EBC {

class PairHMMSampler : public IOptimizable{

	struct is_near {
		bool operator() (pair<double, HMMPathSample>& a, pair<double, HMMPathSample>& b)
		{ return (fabs(b.first - a.first) < 0.000001); }
	};

protected:
	Maths* maths;

	double divergenceT;
	vector<SequenceElement*>* seq1;
	vector<SequenceElement*>* seq2;

	SubstitutionModelBase* substModel;
	IndelModel* indelModel;

	//FIXME - generatie initial model params if you want to use it !
	OptimizedModelParameters modelParams;
	Optimizer* bfgs;

	ForwardPairHMM fwdHmm;
	//ViterbiPairHMM vitHmm;
	BackwardPairHMM bacHmm;

	double forwardLnL;

	unsigned int sampleCount;
	unsigned int sampleMinCount;
	unsigned int sampleMaxCount;
	double sampleDeltaLnL;
	double totalSampleLnL;

	HMMPathSample vitSmpl;

	//log likelihood difference between the best and the worst sample;
	double sampleDelta;

	list<pair<double, HMMPathSample> >samples;

	void sampleInitialSet();

public:
	PairHMMSampler(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, SubstitutionModelBase* smdl,
			IndelModel* imdl, double initialDivergence);

	void doExtraStuff();

	virtual ~PairHMMSampler();

	double runIteration();

	double optimiseDivergenceTime();

	double getDivergence()
	{
		return divergenceT;
	}

	void reSample();
};

} /* namespace EBC */

#endif /* SAMPLING_PAIRHMMSAMPLER_HPP_ */
