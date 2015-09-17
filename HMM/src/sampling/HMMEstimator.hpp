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
 * HMMEstimator.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef HMMESTIMATOR_HPP_
#define HMMESTIMATOR_HPP_


#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include "core/HmmException.hpp"
#include "core/PMatrixTriple.hpp"
#include "core/IOptimizable.hpp"
#include "core/Optimizer.hpp"

#include "models/SubstitutionModelBase.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/GTRModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"

#include "heuristics/GuideTree.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletAligner.hpp"

#include "hmm/ViterbiPairHMM.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "hmm/BackwardPairHMM.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "hmm/BackwardPairHMM.hpp"
#include "hmm/DpMatrixFull.hpp"

#include "sampling/PairHMMSampler.hpp"

#include <sstream>
#include <vector>
#include <array>

using namespace std;


namespace EBC
{

class HMMEstimator : public IOptimizable
{
protected:


	Dictionary* dict;

	SubstitutionModelBase* substModel;
	IndelModel* indelModel;

	Sequences* inputSequences;
	Maths* maths;

	//GuideTree* gtree;
	//TripletSamplingTree tst;

	//TripletAligner* tal;
	ViterbiPairHMM* vphmm;

	OptimizedModelParameters* modelParams;
	Optimizer* bfgs;

	vector<array<unsigned int, 3> > tripletIdxs;

	vector<PairHMMSampler> sampleWorkers;

	unsigned int gammaRateCategories;

	unsigned int tripletIdxsSize;

	bool estimateAlpha;

	double userAlpha;

	void calculateInitialPairs(Definitions::ModelType model, vector<double> substP, vector<double> indelP, double dist);

	double runIteration();

	void optimise();

public:
	HMMEstimator(Sequences* inputSeqs, Definitions::ModelType model,
			Definitions::OptimizationType ot,
			unsigned int rateCategories, double alpha, bool estimateAlpha,
			vector<double> substP, vector<double> indelP, double dist) ;

	virtual ~HMMEstimator();


	vector<double> getSubstitutionParameters();
	vector<double> getIndelParameters();

	double getAlpha();
/*
	GuideTree* getGuideTree()
	{
		return gtree;
	}
*/
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
