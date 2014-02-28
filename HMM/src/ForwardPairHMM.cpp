/*
 * ForwardPairHMM.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "ForwardPairHMM.hpp"

namespace EBC
{

ForwardPairHMM::ForwardPairHMM(Sequences* inputSeqs) :
		EvolutionaryPairHMM(inputSeqs)
{
	initializeModels();
	//initial model calculation
	calculateModels();
	getSequencePair();
	initializeStates();
	setTransitionProbabilities();

}

ForwardPairHMM::~ForwardPairHMM()
{
	// TODO Auto-generated destructor stub
}

void ForwardPairHMM::initializeModels()
{
	generateInitialParameters();
	//start time is the first parameter
	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
}

void ForwardPairHMM::runForwardAlgorithm()
{
	unsigned int i;
	unsigned int j;

	double xx,xy,xm,yx,yy,ym,mx,my,mm;

	double emissionM;
	double emissionX;
	double emissionY;

	while (i != xSize && j != ySize)

	for (i = 1; i<xSize; i++)
	{
		for (j = 1; j<ySize; j++)
		{
			emissionM = log(substModel->getPXiYi(seq1[i].getMatrixIndex(), seq2[j].getMatrixIndex()));
			emissionX = log(substModel->getQXi(seq1[i].getMatrixIndex()));
			emissionY = log(substModel->getQXi(seq2[j].getMatrixIndex()));

			//Mij = PXiYj((Tmm*Mi-1,j-1)+(Txm*Xi-1,j-1)+(Tym*Yi-1,j-1));
			//in log space

			xm = (*M)(i-1,j) + X->getTransitionProbability(M);
			xx = (*X)(i-1,j) + X->getTransitionProbability(X);
			xy = (*Y)(i-1,j) + X->getTransitionProbability(Y);

			ym = (*M)(i,j-1) + Y->getTransitionProbability(M);
			yx = (*X)(i,j-1) + Y->getTransitionProbability(X);
			yy = (*Y)(i,j-1) + Y->getTransitionProbability(Y);

			mm = (*M)(i-1,j-1) + M->getTransitionProbability(M);
			mx = (*X)(i-1,j-1) + M->getTransitionProbability(X);
			my = (*Y)(i-1,j-1) + M->getTransitionProbability(Y);

		}
	}
}

} /* namespace EBC */
