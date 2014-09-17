/*
 * EvolutionaryPairHMM.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#include "hmm/EvolutionaryPairHMM.hpp"
#include "models/NegativeBinomialGapModel.hpp"

namespace EBC
{

EvolutionaryPairHMM::EvolutionaryPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dct,
		unsigned int rateCategories, Maths* mths, Definitions::ModelType model, bool banding, unsigned int bandPercentage, Definitions::DpMatrixType mt)
		: maths(mths), gammaRateCategories(rateCategories)
{
	M = X = Y = NULL;
	mlParameters = NULL;
	this->dict = dct;

	this->seq1 = s1;
	this->seq2 = s2;

	this->xSize = seq1.size() +1;
	this->ySize = seq2.size() +1;

	maths  = new Maths();
	DEBUG("Creating the gap model");
	indelModel = new NegativeBinomialGapModel();

	DEBUG("Creating the model");
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

	bandFactor = bandPercentage;
	bandingEnabled = banding;

	//initialize parameter arrays
	//TODO - set parameters depending on the values provided
	getBandWidth();
	initializeStates(mt);

}

void EvolutionaryPairHMM::setModelParameters(std::vector<double> indel_params,
		std::vector<double> subst_params,double evolDistance, double alpha)
{
	this->indelParameters = indelModel->getParamsNumber();
	this->substParameters = substModel->getParamsNumber();
	this->totalParameters = indelParameters + substParameters +1; //plus time

	if (indelParameters == indel_params.size() && substParameters == subst_params.size())
	{
		indelModel->setParameters(indel_params);
		indelModel->setTime(evolDistance);
		substModel->setAlpha(alpha);
		substModel->setParameters(subst_params);
		substModel->setTime(evolDistance);
	}
	else
	{
		throw HmmException("setModelParameters : wrong number of parameters provided");
	}
}

void EvolutionaryPairHMM::setModelFrequencies(double* freqs)
{
	substModel->setObservedFrequencies(freqs);
}

void EvolutionaryPairHMM::setTransitionProbabilities()
{
	double e,g;

	e = indelModel->getGapExtensionProbability();
	g = indelModel->getGapOpeningProbability();

	M->setTransitionProbabilityFromMatch(log(1-2*g));
	M->setTransitionProbabilityFromInsert(log((1-e)*(1-2*g)));
	M->setTransitionProbabilityFromDelete(log((1-e)*(1-2*g)));

	X->setTransitionProbabilityFromInsert(log(e+((1-e)*g)));
	Y->setTransitionProbabilityFromDelete(log(e+((1-e)*g)));

	X->setTransitionProbabilityFromDelete(log((1-e)*g));
	Y->setTransitionProbabilityFromInsert(log((1-e)*g));

	X->setTransitionProbabilityFromMatch(log(g));
	Y->setTransitionProbabilityFromMatch(log(g));
}

double* EvolutionaryPairHMM::generateInitialSubstitutionParameters()
{
	//don't include time
	double* params = new double[this->substParameters-1];
	for(unsigned int i=0; i< this->substParameters-1; i++)
	{
		params[i] = 0.2 + 0.1*maths->rndu();
	}
	return params;
}
double* EvolutionaryPairHMM::generateInitialIndelParameters()
{
	//don't include time
	double* params = new double[this->indelParameters-1];
	for(unsigned int i=0; i< this->indelParameters-1; i++)
	{
		params[i] = 0.2 + 0.1*maths->rndu();
	}
	return params;
}

double EvolutionaryPairHMM::generateInitialDistanceParameter()
{
	return 0.2 + 0.1*maths->rndu();
}


void EvolutionaryPairHMM::summarize()
{
	double e,g;
	e = indelModel->getGapExtensionProbability();
	g = indelModel->getGapOpeningProbability();

	cout << " Transition probabilities: " << endl;
	cout << "M->M : " << 1-2*g << endl;
	cout << "I->I : " << e+((1-e)*g) << endl;
	cout << "M->I : " << g << endl;
	cout << "I->M : " << (1-2*g)*(1-e) << endl;
	cout << "I->D : " << (1-e)*g << endl << endl;

	indelModel->summarize();
	substModel->summarize();

}

void EvolutionaryPairHMM::initializeStates(Definitions::DpMatrixType mt)
{

	if (M != NULL)
		delete M;
	if (X != NULL)
		delete X;
	if (Y != NULL)
		delete Y;

	switch (mt)
	{
	case Definitions::DpMatrixType::Full :
		M = new PairwiseHmmMatchState(xSize,ySize);
		X = new PairwiseHmmInsertState(xSize,ySize);
		Y = new PairwiseHmmDeleteState(xSize,ySize);
		break;
	case Definitions::DpMatrixType::Limited :
		M = new PairwiseHmmMatchState(new DpMatrixLoMem(xSize,ySize));
		X = new PairwiseHmmInsertState(new DpMatrixLoMem(xSize,ySize));
		Y = new PairwiseHmmDeleteState(new DpMatrixLoMem(xSize,ySize));
		break;
	default :
		M = new PairwiseHmmMatchState(xSize,ySize);
		X = new PairwiseHmmInsertState(xSize,ySize);
		Y = new PairwiseHmmDeleteState(xSize,ySize);
	}
}

void EvolutionaryPairHMM::calculateModels()
{
	substModel->calculateModel();
	indelModel->calculate();
}

EvolutionaryPairHMM::~EvolutionaryPairHMM()
{
	// TODO Auto-generated destructor stub
	//delete bfgs;
	//delete Y;
	//delete X;
	//delete M;
	//delete[] optParameters;
	delete[] mlParameters;
	//delete substModel;
	//delete indelModel;
    delete maths;
}

} /* namespace EBC */
