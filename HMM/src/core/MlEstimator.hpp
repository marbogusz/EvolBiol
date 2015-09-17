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
 * MlEstimator.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef MLESTIMATOR_HPP_
#define MLESTIMATOR_HPP_


#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/MlEstimator.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "hmm/ViterbiPairHMM.hpp"
#include "models/SubstitutionModelBase.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "models/IndelModel.hpp"
#include "core/Sequences.hpp"
#include <sstream>
#include "core/HmmException.hpp"
#include "core/PMatrixDouble.hpp"
#include "core/IOptimizable.hpp"
#include "core/BrentOptimizer.hpp"

#include <vector>
#include <map>
#include <array>

using namespace std;


namespace EBC
{

class MlEstimator : public IOptimizable
{

private:

protected:

	BrentOptimizer* numopt;

	Dictionary* dict;
	SubstitutionModelBase* substModel;
	IndelModel* indelModel;
	Sequences* inputSequences;
	Maths* maths;

	unsigned int gammaRateCategories;

	bool estimateSubstitutionParams;
	bool estimateIndelParams;
	bool estimateDivergence;
	bool estimateAlpha;

	unsigned int pairCount;

	//for viterbi calculation
	vector<ViterbiPairHMM*> hmms;

	unsigned int currentPair;

	vector<double> optTimes;

	//
	//vector<SubstitutionModelBase*> substs;

	vector<PMatrixDouble*> ptMatrices;

	OptimizedModelParameters* modelParams;

	vector<map<array<short, 2>, unsigned int> > patterns;

	bool useViterbi;

public:
	MlEstimator(Sequences* inputSeqs, Definitions::ModelType model,std::vector<double> indel_params,
			std::vector<double> subst_params,Definitions::OptimizationType ot,
			unsigned int rateCategories, double alpha, bool estimateAlpha, vector<double> userTimes, bool useViterbi);

	virtual ~MlEstimator();

	void runIteration(const column_vector& bfgsParameters);

	double runIteration();

	const vector<double> getOptimizedTimes()
	{
		//return this->modelParams->getDivergenceTimes();
		return optTimes;
	}

	//ModelParameters getMlParameters()
	//{
	//	return this->modelParameters;
	//}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
