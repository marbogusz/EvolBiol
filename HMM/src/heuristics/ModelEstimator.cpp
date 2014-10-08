/*
 * ModelEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "heuristics/ModelEstimator.hpp"

namespace EBC
{

ModelEstimator::ModelEstimator(Sequences* inputSeqs, Definitions::ModelType model ,
		Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, bool estimateAlpha) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories),
				gtree(inputSeqs), tst(gtree)
{
	DEBUG("About to sample some triplets");

	vector<array<unsigned int, 3> > tripletIdxs = tst.sampleFromTree();

	for (auto idx : tripletIdxs)
	{
		TripletAligner tal(inputSequences, idx, gtree.getDistanceMatrix());
		tripleAlignments.push_back(tal.align());
	}

	sme = new SubstitutionModelEstimator(inputSeqs, model ,ot, rateCategories, alpha, estimateAlpha, tripletIdxs.size());

	for(int al = 0; al < tripleAlignments.size(); al++)
	{
		sme->addTriplet(tripleAlignments[al]);
	}

	sme->optimize();

	double tb1,tb2,tb3;

	ste = new StateTransitionEstimator(ot);

	for(int al = 0; al < tripleAlignments.size(); al++)
	{
		tb1 = sme->getModelParams()->getDivergenceTime(al*3);
		tb2 = sme->getModelParams()->getDivergenceTime((al*3)+1);
		tb3 = sme->getModelParams()->getDivergenceTime((al*3)+2);
		ste->addPair(tripleAlignments[al][0],tripleAlignments[al][1],tb1+tb2);
		ste->addPair(tripleAlignments[al][1],tripleAlignments[al][2],tb2+tb3);
	}
	ste->optimize();
}

ModelEstimator::~ModelEstimator()
{
    delete maths;
    delete sme;
    delete ste;
}

vector<double> ModelEstimator::getSubstitutionParameters()
{
	return this->sme->getModelParams()->getSubstParameters();
}

vector<double> ModelEstimator::getIndelParameters()
{
	return this->ste->getModelParams()->getIndelParameters();
}

double ModelEstimator::getAlpha()
{
	return this->sme->getModelParams()->getAlpha();
}

} /* namespace EBC */
