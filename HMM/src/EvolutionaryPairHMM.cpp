/*
 * EvolutionaryPairHMM.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#include "EvolutionaryPairHMM.hpp"
#include "GTRModel.hpp"
#include "AffineGeometricGapModel.hpp"
#include "PairHmmMatchState.hpp"
#include "PairHmmInsertionState.hpp"
#include "PairHmmDeletionState.hpp"

namespace EBC
{

EvolutionaryPairHMM::EvolutionaryPairHMM(Sequences* inputSeqs) : inputSequences(inputSeqs)
{
	M = X = Y = NULL;
	dict = inputSeqs->getDictionary();
	maths  = new Maths();
	DEBUG("Creating the substitution model");
	substModel = new GTRModel(dict, maths);
	DEBUG("Creating the gap model");
	indelModel = new AffineGeometricGapModel();
}

void EvolutionaryPairHMM::runForwardAlgorithm()
{

}

void EvolutionaryPairHMM::setTransitionProbabilities()
{
	double e,g;

	e = indelModel->getGapExtensionProbability();
	g = indelModel->getGapOpeningProbability();

	M->setTransitionProbability(M,log(1-2*g));
	M->setTransitionProbability(X,log((1-e)*(1-2*g)));
	M->setTransitionProbability(Y,log((1-e)*(1-2*g)));

	X->setTransitionProbability(X,log(e));
	Y->setTransitionProbability(Y,log(e));

	X->setTransitionProbability(Y,log((1-e)*g));
	Y->setTransitionProbability(X,log((1-e)*g));

	X->setTransitionProbability(M,log(g));
	Y->setTransitionProbability(M,log(g));
}

void EvolutionaryPairHMM::initializeStates()
{
	double e,g;
	e = indelModel->getGapExtensionProbability();
	g = indelModel->getGapOpeningProbability();

	if (M != NULL)
		delete M;
	if (X != NULL)
		delete X;
	if (Y != NULL)
		delete Y;

	M = new PairHmmMatchState(xSize,ySize,e,g);
	X = new PairHmmInsertionState(xSize,ySize,e,g);
	Y = new PairHmmDeletionState(xSize,ySize,e,g);

	M->addTransitionProbability(M,log(1.0-2.0*g));
	M->addTransitionProbability(X,log((1.0-e)*(1.0-2.0*g)));
	M->addTransitionProbability(Y,log((1.0-e)*(1.0-2.0*g)));

	X->addTransitionProbability(X,log(e));
	Y->addTransitionProbability(Y,log(e));

	X->addTransitionProbability(Y,log((1.0-e)*g));
	Y->addTransitionProbability(X,log((1.0-e)*g));

	X->addTransitionProbability(M,log(g));
	Y->addTransitionProbability(M,log(g));
}


void EvolutionaryPairHMM::getSequencePair()
{
	this->seq1 = inputSequences->getSequencesAt(0);
	this->seq2 = inputSequences->getSequencesAt(1);
	this->xSize = seq1.size() +1;
	this->ySize = seq2.size() +1;
}

void EvolutionaryPairHMM::generateInitialParameters()
{
	 //time is a parameter with both indel and subst, we use 1 common time
	this->indelParameters = indelModel->getParamsNumber();
	this->substParameters = substModel->getParamsNumber();
	this->totalParameters = indelParameters + substParameters -1;
	this->mlParameters = new double[totalParameters];

	mlParameters[0] = 0.0; // first parameter hack
	unsigned int sp = this->substParameters-2;
	double tempVal;

	for(unsigned i=1; i< totalParameters; i++)
	{
		tempVal = 0.2 + 0.1*maths->rndu();
		mlParameters[i] = tempVal;//sp>0 ? log(tempVal) : Maths::logit(tempVal);
		sp --;
	}
}


void EvolutionaryPairHMM::calculateModels()
{
	indelModel->setParameters(mlParameters+substParameters-1);
	substModel->setParameters(mlParameters);
	substModel->setParametersInMatrix();
	substModel->setDiagMeans();
	substModel->doEigenDecomposition();
	substModel->calculatePt();
}

EvolutionaryPairHMM::~EvolutionaryPairHMM()
{
	// TODO Auto-generated destructor stub
}

} /* namespace EBC */
