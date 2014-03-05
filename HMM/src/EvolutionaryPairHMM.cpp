/*
 * EvolutionaryPairHMM.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#include "EvolutionaryPairHMM.hpp"
#include "GTRModel.hpp"
#include "HKY85Model.hpp"
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
	substModel = new HKY85Model(dict, maths);
	//substModel = new HKY85Model(dict, maths);
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


	DEBUG ("Opening probs " << g);
	DEBUG ("Extension probs " << e);




	if (M != NULL)
		delete M;
	if (X != NULL)
		delete X;
	if (Y != NULL)
		delete Y;

	M = new PairHmmMatchState(xSize,ySize,g,e);
	X = new PairHmmInsertionState(xSize,ySize,g,e);
	Y = new PairHmmDeletionState(xSize,ySize,g,e);

	M->addTransitionProbabilityFrom(M,log(1.0-2.0*g));
	M->addTransitionProbabilityFrom(X,log((1.0-e)*(1.0-2.0*g)));
	M->addTransitionProbabilityFrom(Y,log((1.0-e)*(1.0-2.0*g)));

	X->addTransitionProbabilityFrom(X,log(e));
	Y->addTransitionProbabilityFrom(Y,log(e));

	X->addTransitionProbabilityFrom(Y,log((1.0-e)*g));
	Y->addTransitionProbabilityFrom(X,log((1.0-e)*g));

	X->addTransitionProbabilityFrom(M,log(g));
	Y->addTransitionProbabilityFrom(M,log(g));
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

	mlParameters[0] = 1; // first parameter hack
	double tempVal;
	for(unsigned i=1; i< totalParameters; i++)
	{
		tempVal = 0.2 + 0.1*maths->rndu();
		mlParameters[i] = tempVal;
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
