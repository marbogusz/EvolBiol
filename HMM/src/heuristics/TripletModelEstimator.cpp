/*
 * TripletModelEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "heuristics/TripletModelEstimator.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletAligner.hpp"
#include "heuristics/StateTransitionEstimator.hpp"

namespace EBC
{

TripletModelEstimator::TripletModelEstimator(Sequences* inputSeqs, Definitions::ModelType model ,
		Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, bool estimateAlpha) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories),
				gtree(inputSeqs), tst(gtree)
{
	maths = new Maths();
	dict = inputSequences->getDictionary();
	//we're running triplets, sequence count is 3!

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



	DEBUG("Gamma rates " << gammaRateCategories);

	estimateSubstitutionParams = true;
	this->estimateAlpha = estimateAlpha;

	DEBUG("About to sample some triplets");
	vector<array<unsigned int, 3> > tripletIdxs = tst.sampleFromTree();

	patterns.resize(tripletIdxs.size());
	ptMatrices.resize(tripletIdxs.size());

	DEBUG("Creating branch-specific substitution models");

	for(unsigned int i =0; i<ptMatrices.size(); i++)
	{
		ptMatrices[i][0] = new PMatrixTriple(substModel);
		ptMatrices[i][1] = new PMatrixTriple(substModel);
		ptMatrices[i][2] = new PMatrixTriple(substModel);
	}

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies(tripletIdxs[0]));

	for (auto idx : tripletIdxs)
	{
		TripletAligner tal(inputSequences, idx, gtree.getDistanceMatrix());
		tripleAlignments.push_back(tal.align());
	}

	for(int al = 0; al < tripleAlignments.size(); al++)
	{
		for(int pos = 0; pos < tripleAlignments[al][0].size(); pos++)
		{
			//count patterns
			patterns[al][{{tripleAlignments[al][0][pos].getMatrixIndex(), tripleAlignments[al][1][pos].getMatrixIndex(),tripleAlignments[al][2][pos].getMatrixIndex()}}]++;
		}
	}

	//TripletSamplingTree tst(GuideTree(inputSequences));
	//triplets = tst.sampleFromDM();

	modelParams = new OptimizedModelParameters(substModel, NULL,3, 3*tripleAlignments.size(), estimateSubstitutionParams,
			false, estimateAlpha, true, maths);

	modelParams->setAlpha(alpha);

	bfgs = new Optimizer(modelParams, this,ot);
	bfgs->optimize();
	//we now have the triplet times + alignments
	//optimize indels - construct a vector of alignment pairs with the distance

	modelParams->outputParameters();

	StateTransitionEstimator* ste = new StateTransitionEstimator(ot);

	double tb1,tb2,tb3;

	for(int al = 0; al < tripleAlignments.size(); al++)
	{
		tb1 = modelParams->getDivergenceTime(al*3);
		tb2 = modelParams->getDivergenceTime((al*3)+1);
		tb3 = modelParams->getDivergenceTime((al*3)+2);
		ste->addPair(tripleAlignments[al][0],tripleAlignments[al][1],tb1+tb2);
		ste->addPair(tripleAlignments[al][1],tripleAlignments[al][2],tb2+tb3);
	}

	ste->optimize();

}

TripletModelEstimator::~TripletModelEstimator()
{
	// TODO Auto-generated destructor stub
	delete bfgs;
	delete modelParams;
    delete maths;
}

vector<double> TripletModelEstimator::getSubstitutionParameters()
{
	return this->modelParams->getSubstParameters();
}

vector<double> TripletModelEstimator::getIndelParameters()
{
	return this->modelParams->getIndelParameters();
}

double TripletModelEstimator::getAlpha()
{
	return this->modelParams->getAlpha();
}

double TripletModelEstimator::runIteration()
{
	double result = 0;


	double partial1;

	substModel->setAlpha(modelParams->getAlpha());
	substModel->setParameters(modelParams->getSubstParameters());
	substModel->calculateModel();


	//substModel->summarize();

	for (unsigned int i = 0; i< ptMatrices.size(); i++)
	{
		for(unsigned int j=0;j<3;j++)
		{
			ptMatrices[i][j]->setTime(modelParams->getDivergenceTime(3*i +j));
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
			}  //->getPairSitePattern(it.first[0],it.first[1]) * it.second;
			result += log(partial1)* it.second;
			//cout << "Pattern " << it.first[0] << " "  << it.first[1] << " " << it.first[2]  << "\t\t\t" << it.second << "\t\t" << log(partial1) << endl;
		}
	}
	//DEBUG("lnl result:" << result);
	return result * -1.0;

}

} /* namespace EBC */
