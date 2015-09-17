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
 * MlEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "core/MlEstimator.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"

namespace EBC
{

MlEstimator::MlEstimator(Sequences* inputSeqs, Definitions::ModelType model ,std::vector<double> indel_params,
		std::vector<double> subst_params,   Definitions::OptimizationType ot,
		unsigned int rateCategories, double alpha, bool estimateAlpha, vector<double> userTime, bool alignViterbi) : inputSequences(inputSeqs), gammaRateCategories(rateCategories),
		pairCount(inputSequences->getPairCount()), hmms(pairCount), ptMatrices(pairCount, nullptr), useViterbi(alignViterbi), patterns(pairCount),
		optTimes(pairCount)

{
	maths = new Maths();
	dict = inputSequences->getDictionary();

	if (model == Definitions::ModelType::GTR)
	{
		substModel = new GTRModel(dict, maths,gammaRateCategories);
	}
	else if (model == Definitions::ModelType::HKY85)
	{
		substModel = new HKY85Model(dict, maths,gammaRateCategories);
	}
	else if (model == Definitions::ModelType::LG)
	{
			substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaLgModel);
	}

	indelModel = new NegativeBinomialGapModel();

	substModel->setAlpha(alpha);

	/*
	estimateSubstitutionParams = true;
	if (subst_params.size() > 0){
		estimateSubstitutionParams = false;
	}
	*/

	estimateIndelParams = false;
	estimateSubstitutionParams = false;
	//can't estimate alpha!
	//unless done before using triplets!

	this->estimateAlpha = false;

	//hack here!

	useViterbi = false;

	/*
	modelParams = new OptimizedModelParameters(substModel, NULL,inputSequences->getSequenceCount(), pairCount, estimateSubstitutionParams,
			estimateIndelParams, estimateAlpha, true, maths);

	*/
	modelParams = new OptimizedModelParameters(substModel, NULL,2, 1, false,
				false, false, true, maths);

	if (estimateSubstitutionParams)
		modelParams->generateInitialSubstitutionParameters();
	//modelParams->generateInitialDistanceParameters();

	//modelParams->setUserDivergenceParams(userTime);

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
	if (!estimateSubstitutionParams){
		substModel->setParameters(subst_params);
		substModel->calculateModel();
	}

	//bfgs = new BFGS(this,ot);
	numopt = new BrentOptimizer(modelParams, NULL);

	for(unsigned int i =0; i<pairCount; i++)
	{
		ptMatrices[i] = new PMatrixDouble(substModel);
		std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);
		vector<SequenceElement*>*  s1 = inputSequences->getSequencesAt(idxs.first);
		vector<SequenceElement*>*  s2 = inputSequences->getSequencesAt(idxs.second);
		for(int j = 0; j< s1->size(); j++)
		{
			patterns[i][{{(*s1)[j]->getMatrixIndex(),(*s2)[j]->getMatrixIndex()}}]++;
		}
		currentPair = i;
		modelParams->setUserDivergenceParams({userTime[i]});
		INFO("Running pairwise calculator for sequence id " << idxs.first << " and " << idxs.second);

		numopt->setTarget(this);
		numopt->setAccuracy(1e-6);
		numopt->setBounds(0.0000001, modelParams->divergenceBound);

		numopt->optimize() * -1.0;

		//bfgs->optimize();
		optTimes[i] = modelParams->getDivergenceTime(0);

	}
}

MlEstimator::~MlEstimator()
{
	//delete bfgs;
	delete numopt;
	delete modelParams;
    delete maths;
}

double MlEstimator::runIteration()
{
	double result = 0;
	if (estimateSubstitutionParams){
		substModel->setParameters(modelParams->getSubstParameters());
		substModel->calculateModel();
	}
	//for(unsigned int i =0; i<pairCount; i++)
	//{
		//this calculates the matrix(matrices for a gamma model)
		ptMatrices[currentPair]->setTime(modelParams->getDivergenceTime(0));
		ptMatrices[currentPair]->calculate();

		//go through the map of patterns!
		for(auto it : patterns[currentPair])
		{
			result += ptMatrices[currentPair]->getPairSitePattern(it.first[0],it.first[1]) * it.second;
		}

	//}

	//cerr << result << endl;
	return result * -1.0;
}


} /* namespace EBC */
