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
 * BandingEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "core/BandingEstimator.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"
#include "hmm/DpMatrixFull.hpp"
#include "extras/LikelihoodSurfacePlotter.hpp"

namespace EBC
{

BandingEstimator::BandingEstimator(Definitions::AlgorithmType at, Sequences* inputSeqs, Definitions::ModelType model ,std::vector<double> indel_params,
		std::vector<double> subst_params, Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, GuideTree* g) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories), pairCount(inputSequences->getPairCount()),
				/*hmms(pairCount), bands(pairCount),*/ divergenceTimes(pairCount), algorithm(at)//, gt(g)
{
	//Banding estimator means banding enabled!

	DEBUG("Starting Banding Estimator");
	maths = new Maths();
	dict = inputSequences->getDictionary();

	//Helper models
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

	estimateSubstitutionParams = false;
	estimateIndelParams = false;
	estimateAlpha = false;
	/*
	estimateSubstitutionParams = subst_params.size() != substModel->getParamsNumber();
	estimateIndelParams = indel_params.size() == 0;
	//Do not estimate alpha here. Alpha needs to be estimated beforehand
	this->estimateAlpha = false;

	DEBUG("Pairwise banding model estimator starting");
	DEBUG("Estimate substitution parameters set to : " << estimateSubstitutionParams << " Estimate indel parameters set to : " << estimateIndelParams);
	DEBUG("Estimate alpha set to : " << estimateAlpha << " , rate categories " << gammaRateCategories << " , alpha : " << alpha);

	FileLogger::DebugLogger() << "Estimate substitution parameters set to : " << estimateSubstitutionParams << " Estimate indel parameters set to : " << estimateIndelParams << "\n";
	FileLogger::DebugLogger() << "Estimate alpha set to : " << estimateAlpha << " , rate categories " << gammaRateCategories << " , alpha : " << alpha << "\n";
	 */
	//pairwise comparison mode
	modelParams = new OptimizedModelParameters(substModel, indelModel,2, 1, estimateSubstitutionParams,
			estimateIndelParams, estimateAlpha, true, maths);

	modelParams->boundDivergenceBasedOnLambda(indel_params[0]);

	if(!estimateIndelParams)
		modelParams->setUserIndelParams(indel_params);
	if(!estimateSubstitutionParams)
		modelParams->setUserSubstParams(subst_params);
	modelParams->setAlpha(alpha);

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
	if (estimateSubstitutionParams == false)
	{
		//set parameters and calculate the model
		substModel->setAlpha(modelParams->getAlpha());
		substModel->setParameters(modelParams->getSubstParameters());
		substModel->calculateModel();
	}

	if (estimateIndelParams == false)
	{
			//set parameters and calculate the model
		indelModel->setParameters(modelParams->getIndelParameters());
	}

	//let's assume that we have all the parameters estimated
	//need to get times!

	EvolutionaryPairHMM *hmm;
//non banden probs
//	vector<EvolutionaryPairHMM*> hmmsNB(pairCount);

/*
	for(unsigned int i =0; i< pairCount; i++)
		{
			std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);
			DEBUG("Running band calculator for sequence " << idxs.first << " and " << idxs.second);
			FileLogger::InfoLogger() << "Running band calculator for sequence " << idxs.first << " and " << idxs.second << "\n";
			BandCalculator* bc = new BandCalculator(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
					substModel, indelModel, gt->getDistanceMatrix()->getDistance(idxs.first,idxs.second));
			bands[i] = bc->getBand();
			if (at == Definitions::AlgorithmType::Viterbi)
			{
				hmm = hmms[i] = new ViterbiPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
					substModel, indelModel, Definitions::DpMatrixType::Full, bands[i]);
			}
			else if (at == Definitions::AlgorithmType::Forward)
			{
				hmm = hmms[i] = new ForwardPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
						substModel, indelModel, Definitions::DpMatrixType::Full, bands[i]);
				//hmmsNB[i] = new ForwardPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
				//						substModel, indelModel, Definitions::DpMatrixType::Full, NULL);
			}
			else
			{
				throw HmmException("Wrong algorithm type - use either Forward or viterbi\n");
			}
			delete bc;
		}
*/
	numopt = new BrentOptimizer(modelParams, NULL);
	//bfgs->optimize();
	//this->modelParams->logParameters();

}

BandingEstimator::~BandingEstimator()
{
	//for(auto hmm : hmms)
	//	delete hmm;
	delete numopt;
	delete modelParams;
    delete maths;
    //for (auto bnd : bands)
    //	delete bnd;
    delete indelModel;
    delete substModel;
}

void BandingEstimator::optimizePairByPair()
{
	ForwardPairHMM* fhmm;
	BackwardPairHMM* bhmm;
	Band* band;
	//DistanceMatrix* dm = gt->getDistanceMatrix();
	PairHmmCalculationWrapper* wrapper = new PairHmmCalculationWrapper();
	double result1, result2;

	for(unsigned int i =0; i< pairCount; i++)
	//int i = 16;
	{
		DEBUG("Optimizing distance for pair #" << i);
		std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);
		INFO("Running pairwise calculator for sequence id " << idxs.first << " and " << idxs.second
				<< " ,number " << i+1 <<" out of " << pairCount << " pairs" );
		//BandCalculator* bc = new BandCalculator(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
		//		substModel, indelModel, gt->getDistanceMatrix()->getDistance(idxs.first,idxs.second));
		band = NULL;//bc->getBand();

		DEBUG("Creating forward algorithm to optimize the pairwise divergence time...");
		fhmm = new ForwardPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
				substModel, indelModel, Definitions::DpMatrixType::Full, band);
		bhmm = new BackwardPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
				substModel, indelModel, Definitions::DpMatrixType::Full, band);


		//hmm->setDivergenceTimeAndCalculateModels(modelParams->getDivergenceTime(0)); //zero as there's only one pair!

		//LikelihoodSurfacePlotter lsp;
		//lsp.setTargetHMM(hmm);
		//lsp.getLikelihoodSurface();

/*
		wrapper->setTargetHMM(fhmm);
		DUMP("Set model parameter in the hmm...");
		wrapper->setModelParameters(modelParams);
		modelParams->setUserDivergenceParams({0.2});
		numopt->setTarget(wrapper);
		numopt->setAccuracy(Definitions::highDivergenceAccuracyDelta);
		//FIXME - hardcoded right bound
		numopt->setBounds(Definitions::almostZero, modelParams->divergenceBound);


		result = numopt->optimize() * -1.0;
*/
		DEBUG("*#*#*#*#*# Final forward and backward calculation for divergence " << modelParams->getDivergenceTime(0));

		fhmm->setDivergenceTimeAndCalculateModels(0.5);//modelParams->getDivergenceTime(0));
		result1 = fhmm->runAlgorithm();
		bhmm->setDivergenceTimeAndCalculateModels(0.5);//modelParams->getDivergenceTime(0));
		result2 = bhmm->runAlgorithm();
		bhmm->calculatePosteriors(fhmm);

		//bhmm->calculateMaximumPosteriorMatrix();

		//auto mp1 = bhmm->getMPAlignment();


		//DUMP("MPD aligment for sequence id " << idxs.first << " and " << idxs.second);
		//DUMP(mp1.first);
		//DUMP(mp1.second);


		cout << inputSequences->getSequenceName(idxs.first) << "\t" <<  inputSequences->getSequenceName(idxs.second) << "\t" << bhmm->getAlignmentPosteriors(inputSequences->getAlignmentsAt(idxs.first), inputSequences->getAlignmentsAt(idxs.second));

		//DEBUG("Likelihood after pairwise optimization: " << result);
		if (result1 <= (Definitions::minMatrixLikelihood /2.0) || result2 <= (Definitions::minMatrixLikelihood /2.0))
		{
			ERROR("Optimization failed for pair #" << i << " Zero probability FWD");
		}
		this->divergenceTimes[i] = modelParams->getDivergenceTime(0);

		if (band != NULL)
			delete band;
		//delete bc;
		delete fhmm;
		delete bhmm;

	}

	DEBUG("Optimized divergence times:");
	DEBUG(this->divergenceTimes);
}


double BandingEstimator::runIteration()
{
	double result = 0;
	double tmp;
	return result;
}

void BandingEstimator::outputDistanceMatrix(stringstream& ss)
{
	unsigned int count, pairCount;
	count = this->inputSequences->getSequenceCount();
	pairCount = this->inputSequences->getPairCount();

	ss << "\t" << this->inputSequences->getSequenceCount() << endl;

	for (unsigned int i = 0; i< count; i++)
	{
		ss << "S" << i << " ";
		for(unsigned int j=0; j< count; j++)
		{
			ss << this->modelParams->getDistanceBetween(i,j) << " ";
		}
		ss << endl;
	}
}

} /* namespace EBC */

