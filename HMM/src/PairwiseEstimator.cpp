/*
 * PairwiseEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "PairwiseEstimator.hpp"

namespace EBC
{


PairwiseEstimator::BFGS::BFGS(PairwiseEstimator* enclosing, Definitions::OptimizationType ot) : optimizationType(ot)
{
	parent = enclosing;
	paramsCount = parent->modelParams.optParamCount();
	this->initParams.set_size(paramsCount);
	this->lowerBounds.set_size(paramsCount);
	this->upperBounds.set_size(paramsCount);

	parent->modelParams.toDlibVector(initParams,lowerBounds,upperBounds);

	DEBUG("DLIB optimizer init with " << paramsCount << " parameters");
}

PairwiseEstimator::BFGS::~BFGS()
{
}

double PairwiseEstimator::BFGS::objectiveFunction(const column_vector& bfgsParameters)
{
	this->parent->modelParams.fromDlibVector(bfgsParameters);
	return parent->runIteration();
}


const column_vector PairwiseEstimator::BFGS::objectiveFunctionDerivative(const column_vector& bfgsParameters)
{
	column_vector results(this->paramsCount);
	return results;
}


void PairwiseEstimator::BFGS::optimize()
{
	using std::placeholders::_1;
	std::function<double(const column_vector&)> f_objective= std::bind( &PairwiseEstimator::BFGS::objectiveFunction, this, _1 );

	switch(optimizationType)
	{
		case Definitions::OptimizationType::BFGS:
		{
			dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
					dlib::objective_delta_stop_strategy(1e-8),
					f_objective,
					derivative(f_objective),
					initParams,
					lowerBounds,
					upperBounds);
			break;
		}
		case Definitions::OptimizationType::BOBYQA:
		{
			dlib::find_min_bobyqa(f_objective, initParams, 10,
					lowerBounds,upperBounds, 0.05, 1e-6, 10000 );
			break;
		}
	}
	DEBUG("BFGS return: " << initParams );
}


PairwiseEstimator::PairwiseEstimator(Sequences* inputSeqs, Definitions::ModelType model ,std::vector<double> indel_params,
		std::vector<double> subst_params, Definitions::OptimizationType ot, bool banding, unsigned int bandPercentage,
		unsigned int rateCategories, double alpha, bool estimateAlpha) : inputSequences(inputSeqs), gammaRateCategories(rateCategories)
{
	this->pairCount = inputSequences->getPairCount();

	this->hmms(pairCount);

	maths = new Maths();
	dict = inputSequences->getDictionary();

	//Helper models
	//FIXME - get some static definitions or sth!!
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

	estimateSubstitutionParams = subst_params.size() == 0;
	estimateIndelParams = indel_params.size() == 0;
	this->estimateAlpha = estimateAlpha;

	this->modelParams(substModel, indelModel, pairCount, estimateSubstitutionParams,
			estimateIndelParams, estimateAlpha, maths);

	if(!estimateIndelParams)
		modelParams.setUserIndelParams(indel_params);
	if(!estimateSubstitutionParams)
		modelParams.setUserSubstParams(subst_params);
	if(!estimateAlpha)
		modelParams.setAlpha(alpha);


	bandFactor = bandPercentage;
	bandingEnabled = banding;

	ForwardPairHMM* hmm;

	for(unsigned int i =0; i<pairCount; i++)
	{
		hmm = hmms[i] = new ForwardPairHMM(inputSequences->getSequencesAt(0), inputSequences->getSequencesAt(1),
				dict, ot , banding, bandPercentage,rateCategories, alpha, maths);
		hmm->setModelFrequencies(inputSequences->getElementFrequencies());
	}

	bfgs = new BFGS(this,ot);
	bfgs->optimize();
}

PairwiseEstimator::~PairwiseEstimator()
{
	// TODO Auto-generated destructor stub
	delete bfgs;
	//delete Y;
	//delete X;
	//delete M;
	//delete substModel;
	//delete indelModel;
    delete maths;
}

double PairwiseEstimator::runIteration()
{
	double result = 0;
	ForwardPairHMM* hmm;

	for(unsigned int i =0; i<pairCount; i++)
	{
		hmm = hmms[i];
		hmm>setModelParameters(modelParams.getIndelParameters(),modelParams.getSubstParameters(),
				modelParams.getDivergenceTime(i), modelParams.getAlpha());
		result += hmm->runForwardAlgorithm();
	}
}

} /* namespace EBC */
