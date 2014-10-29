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

namespace EBC
{

BandingEstimator::BandingEstimator(Definitions::AlgorithmType at, Sequences* inputSeqs, Definitions::ModelType model ,std::vector<double> indel_params,
		std::vector<double> subst_params, Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, GuideTree* gt) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories), pairCount(inputSequences->getPairCount()), hmms(pairCount)
{
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

	estimateSubstitutionParams = subst_params.size() != substModel->getParamsNumber();
	estimateIndelParams = indel_params.size() == 0;
	this->estimateAlpha = estimateAlpha;

	DEBUG("Pairwise banding model estimator starting");
	DEBUG("Estimate substitution parameters set to : " << estimateSubstitutionParams << " Estimate indel parameters set to : " << estimateIndelParams);
	DEBUG("Estimate alpha set to : " << estimateAlpha << " , rate categories " << gammaRateCategories << " , alpha : " << alpha);

	modelParams = new OptimizedModelParameters(substModel, indelModel,inputSequences->getSequenceCount(), pairCount, estimateSubstitutionParams,
			estimateIndelParams, estimateAlpha, true, maths);

	if(!estimateIndelParams)
		modelParams->setUserIndelParams(indel_params);
	if(!estimateSubstitutionParams)
		modelParams->setUserSubstParams(subst_params);
	modelParams->setAlpha(alpha);

	EvolutionaryPairHMM* hmm;

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

	for(unsigned int i =0; i<pairCount; i++)
		{
			std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);
			DEBUG("Running band calculator for sequence " << idxs.first << " and " << idxs.second);
			BandCalculator* bc = new BandCalculator(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
					substModel, indelModel, gt->getDistanceMatrix()->getDistance(idxs.first,idxs.second));
			delete bc;
		}

/*
	for(unsigned int i =0; i<pairCount; i++)
	{
		std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);

		if (at == Definitions::AlgorithmType::Viterbi)
		{
			hmm = hmms[i] = new ViterbiPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
				false, substModel, indelModel, 0, Definitions::DpMatrixType::Limited);
		}
		else if (at == Definitions::AlgorithmType::Forward)
		{
			hmm = hmms[i] = new ForwardPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
					false, substModel, indelModel, 0, Definitions::DpMatrixType::Limited);
		}
		else
		{
			throw HmmException("Wrong algorithm type - use either Forward or viterbi\n");
		}
	}

	bfgs = new Optimizer(modelParams, this, ot);
	bfgs->optimize();
*/
}

BandingEstimator::~BandingEstimator()
{
	// TODO Auto-generated destructor stub
	delete bfgs;
	delete modelParams;
    delete maths;
}

double BandingEstimator::runIteration()
{
	double result = 0;
	EvolutionaryPairHMM* hmm;

	//this->modelParams->outputParameters();
	//cerr << "iteration " << endl;

	if (estimateSubstitutionParams == true)
	{
			//set parameters and calculate the model
		substModel->setAlpha(modelParams->getAlpha());
		substModel->setParameters(modelParams->getSubstParameters());
		substModel->calculateModel();
	}

	if (estimateIndelParams == true)
	{
		//set parameters and calculate the model
		indelModel->setParameters(modelParams->getIndelParameters());
	}



	for(unsigned int i =0; i<pairCount; i++)
	{
		hmm = hmms[i];
		//FIXME - individual indel models!!! or gap opening extension class aggregator!
		// the following is not thread safe for indels!!!
		hmm->setDivergenceTime(modelParams->getDivergenceTime(i));
		indelModel->setTime(modelParams->getDivergenceTime(i));
		indelModel->calculate();
		result += hmm->runAlgorithm();
		//modelParams->outputParameters();
	}
	cerr << " lnl " << result << endl;
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
