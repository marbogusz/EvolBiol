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


MlEstimator::BFGS::BFGS(MlEstimator* enclosing, Definitions::OptimizationType ot) : optimizationType(ot)
{
	parent = enclosing;
	paramsCount = parent->modelParams->optParamCount();
	this->initParams.set_size(paramsCount);
	this->lowerBounds.set_size(paramsCount);
	this->upperBounds.set_size(paramsCount);

	parent->modelParams->toDlibVector(initParams,lowerBounds,upperBounds);

	cerr << "DLIB optimizer init with " << paramsCount << " parameters" << endl;
}

MlEstimator::BFGS::~BFGS()
{
}

double MlEstimator::BFGS::objectiveFunction(const column_vector& bfgsParameters)
{
	this->parent->modelParams->fromDlibVector(bfgsParameters);
	return parent->runIteration();
}


const column_vector MlEstimator::BFGS::objectiveFunctionDerivative(const column_vector& bfgsParameters)
{
	column_vector results(this->paramsCount);
	return results;
}


void MlEstimator::BFGS::optimize()
{
	using std::placeholders::_1;
	std::function<double(const column_vector&)> f_objective= std::bind( &MlEstimator::BFGS::objectiveFunction, this, _1 );
	double likelihood;

	switch(optimizationType)
	{
		case Definitions::OptimizationType::BFGS:
		{
			dlib::objective_delta_stop_strategy strt(1e-7);
			strt.be_verbose();
			likelihood = dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
					strt,
					f_objective,
					derivative(f_objective),
					initParams,
					lowerBounds,
					upperBounds);
			break;
		}
		case Definitions::OptimizationType::BOBYQA:
		{
			likelihood = dlib::find_min_bobyqa(f_objective, initParams, parent->modelParams->optParamCount()+4,
					lowerBounds,upperBounds, 0.05, 1e-7, 20000 );
			break;
		}
	}
	this->parent->modelParams->fromDlibVector(initParams);
	parent->modelParams->logParameters();
}


MlEstimator::MlEstimator(Sequences* inputSeqs, Definitions::ModelType model ,std::vector<double> indel_params,
		std::vector<double> subst_params,   Definitions::OptimizationType ot,
		unsigned int rateCategories, double alpha, bool estimateAlpha, vector<double> userTime, bool alignViterbi) : inputSequences(inputSeqs), gammaRateCategories(rateCategories),
		pairCount(inputSequences->getPairCount()), hmms(pairCount), ptMatrices(pairCount, nullptr), useViterbi(alignViterbi), patterns(pairCount)

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

	estimateSubstitutionParams = true;
	estimateIndelParams = false;
	//can't estimate alpha!
	//unless done before using triplets!

	this->estimateAlpha = false;

	//hack here!

	useViterbi = false;

	modelParams = new OptimizedModelParameters(substModel, NULL,inputSequences->getSequenceCount(), pairCount, true,
				false, false, true, maths);

	modelParams->generateInitialSubstitutionParameters();
	//modelParams->generateInitialDistanceParameters();
	modelParams->setUserDivergenceParams(userTime);

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());

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
	}



	bfgs = new BFGS(this,ot);
	bfgs->optimize();
}

MlEstimator::~MlEstimator()
{
	delete bfgs;
	delete modelParams;
    delete maths;
}

double MlEstimator::runIteration()
{
	double result = 0;

	substModel->setParameters(modelParams->getSubstParameters());
	substModel->calculateModel();

	for(unsigned int i =0; i<pairCount; i++)
	{
		//this calculates the matrix(matrices for a gamma model)
		ptMatrices[i]->setTime(modelParams->getDivergenceTime(i));
		ptMatrices[i]->calculate();

		//go through the map of patterns!
		for(auto it : patterns[i])
		{
			result += ptMatrices[i]->getPairSitePattern(it.first[0],it.first[1]) * it.second;
		}

	}

	//cerr << result << endl;
	return result * -1.0;
}


} /* namespace EBC */
