/*
 * EvolutionaryPairHMM.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#include "EvolutionaryPairHMM.hpp"
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
	DEBUG("Creating the gap model");
	indelModel = new AffineGeometricGapModel();

	//substitution model created by the parent class
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

void EvolutionaryPairHMM::getSequencePair()
{
	this->seq1 = inputSequences->getSequencesAt(0);
	this->seq2 = inputSequences->getSequencesAt(1);
	this->xSize = seq1.size() +1;
	this->ySize = seq2.size() +1;
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
