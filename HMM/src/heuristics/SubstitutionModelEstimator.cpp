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

SubstitutionModelEstimator::SubstitutionModelEstimator(Sequences* inputSeqs, SubstitutionModelBase* model,
		Definitions::OptimizationType ot,unsigned int rateCategories, double alpha,
		bool estimateAlpha, unsigned int matCount) :
				inputSequences(inputSeqs), substModel(model), gammaRateCategories(rateCategories), patterns(matCount), ptMatrices(matCount)


{

	DEBUG("Starting Substitution Model Estimator (SME)");
	DUMP("SME estimate alpha : " << estimateAlpha << " alpha value " << alpha);
	DUMP("SME rate categories : " << rateCategories << " triplets number " << matCount);

	this->estimateSubstitutionParams = true;
	this->estimateAlpha = estimateAlpha;
	this->alpha = alpha;

	currentTriplet = 0;

	maths = new Maths();
	dict = inputSequences->getDictionary();

	modelParams = new OptimizedModelParameters(substModel, NULL,3, 3*patterns.size(), estimateSubstitutionParams,
			false, estimateAlpha, true, maths);

	modelParams->setAlpha(alpha);

	for(int i = 0; i < ptMatrices.size(); i++){
		DUMP("SME: creating ptMatrix");
		ptMatrices[i][0]  = new PMatrixTriple(substModel);
		ptMatrices[i][1]  = new PMatrixTriple(substModel);
		ptMatrices[i][2]  = new PMatrixTriple(substModel);
	}

	bfgs = new Optimizer(modelParams, this,ot);
}

void SubstitutionModelEstimator::clean()
{
	//for(auto entry : ptMatrices)
	//{
	//	delete entry[0];
	//	delete entry[1];
	//	delete entry[2];
	//}


	for (auto itp = patterns.begin(); itp != patterns.end(); itp++)
	{
		for (auto itMap = itp->begin(); itMap != itp->end(); itMap++)
		{
			itMap->second = 0.0;
		}
	}
}

SubstitutionModelEstimator::~SubstitutionModelEstimator()
{
	delete bfgs;
	delete modelParams;
	//delete substModel;
	delete maths;

	for(auto entry : ptMatrices)
	{
		delete entry[0];
		delete entry[1];
		delete entry[2];
	}
}

void SubstitutionModelEstimator::addTriplet(array<vector<unsigned char>*, 3>* tripleAlignment, unsigned int trp, double weight)
{
	DUMP("SME : adding patterns for triplet " << trp << "\t triplet weight : " << weight);
	for(int pos = 0; pos < (*tripleAlignment)[0]->size(); pos++)
	{
		patterns[trp][{{(*(*tripleAlignment)[0])[pos], (*(*tripleAlignment)[1])[pos],(*(*tripleAlignment)[2])[pos]}}] += weight;
	}
	//for (auto pat : patterns[trp])
	//{
	//	DUMP("SME " << pat.first[0] << " " << pat.first[1] << " " << pat.first[2] << " : " << pat.second);
	//}
}

void SubstitutionModelEstimator::optimize()
{

	bfgs->optimize();
	INFO("SubstitutionModelEstimator results:");

	modelParams->logParameters();
/*
	cerr << "kappa\tlnL" << endl;

	double kappa = 0.5;
	double result;
	double partial1;
	while (kappa < 5.0)
	{
		result = 0;
		substModel->setAlpha(modelParams->getAlpha());
		substModel->setParameters({kappa});
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
		cerr << kappa << "\t" << result << endl;
		kappa += 0.025;
	}
*/

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
	//cerr << modelParams->getSubstParameters()[0] << "\t" << result << endl;

	return result * -1.0;



}

} /* namespace EBC */
