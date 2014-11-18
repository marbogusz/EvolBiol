/*
 * SubstitutionModelEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "heuristics/SubstitutionModelEstimator.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletAligner.hpp"
#include "heuristics/StateTransitionEstimator.hpp"

namespace EBC
{

SubstitutionModelEstimator::SubstitutionModelEstimator(Sequences* inputSeqs, Definitions::ModelType model,
		Definitions::OptimizationType ot,unsigned int rateCategories, double alpha,
		bool estimateAlpha, unsigned int sampleCount) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories), patterns(sampleCount), ptMatrices(sampleCount)


{

	DEBUG("Starting Substitution Model Estimator");

	this->estimateSubstitutionParams = true;
	this->estimateAlpha = estimateAlpha;
	this->alpha = alpha;

	currentTriplet = 0;

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

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());

	modelParams = new OptimizedModelParameters(substModel, NULL,3, 3*patterns.size(), estimateSubstitutionParams,
			false, estimateAlpha, true, maths);

	modelParams->setAlpha(alpha);

	bfgs = new Optimizer(modelParams, this,ot);
}

SubstitutionModelEstimator::~SubstitutionModelEstimator()
{
	delete bfgs;
	delete modelParams;
	delete substModel;
	delete maths;

	for(auto entry : ptMatrices)
	{
		delete entry[0];
		delete entry[1];
		delete entry[2];
	}
}

void SubstitutionModelEstimator::addTriplet(array<vector<SequenceElement>, 3>& tripleAlignment)
{
	for(int pos = 0; pos < tripleAlignment[0].size(); pos++)
	{
		patterns[currentTriplet][{{tripleAlignment[0][pos].getMatrixIndex(), tripleAlignment[1][pos].getMatrixIndex(),tripleAlignment[2][pos].getMatrixIndex()}}]++;
	}
	ptMatrices[currentTriplet][0] = new PMatrixTriple(substModel);
	ptMatrices[currentTriplet][1] = new PMatrixTriple(substModel);
	ptMatrices[currentTriplet][2] = new PMatrixTriple(substModel);

	currentTriplet++;
}

void SubstitutionModelEstimator::optimize()
{
	bfgs->optimize();
	INFO("SubstitutionModelEstimator results:");
	modelParams->logParameters();

}

double SubstitutionModelEstimator::runIteration()
{
	double result = 0;

	double partial1;

	substModel->setAlpha(modelParams->getAlpha());
	substModel->setParameters(modelParams->getSubstParameters());
	substModel->calculateModel();

	for (unsigned int i = 0; i< ptMatrices.size(); i++)
	{
		for(unsigned int j=0;j<Definitions::heuristicsTreeSize;j++)
		{
			ptMatrices[i][j]->setTime(modelParams->getDivergenceTime(Definitions::heuristicsTreeSize*i +j));
			ptMatrices[i][j]->calculate();
		}
	}

	for(int al = 0; al < patterns.size(); al++)
	{
		for(auto it : patterns[al])
		{
			partial1 = 0;
			for(int rt = 0; rt < dict->getAlphabetSize(); rt++)
			{
				partial1 += ptMatrices[al][0]->getTripleSitePattern(rt,it.first, ptMatrices[al][1],ptMatrices[al][2]);
			}
			result += log(partial1)* it.second;

		}
	}
	//DEBUG("lnl result:" << result);
	return result * -1.0;

}

} /* namespace EBC */
