/*
 * ViterbiPairHMM.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "ViterbiPairHMM.hpp"

namespace EBC
{

ViterbiPairHMM::ViterbiPairHMM(Sequences* inputSeqs) :
		EvolutionaryPairHMM(inputSeqs)
{
	// TODO Auto-generated constructor stub
	getSequencePair();

}

ViterbiPairHMM::~ViterbiPairHMM()
{
	// TODO Auto-generated destructor stub
}

void ViterbiPairHMM::initializeModels()
{
	generateInitialParameters();
	//start time is the first parameter
	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
}

double ViterbiPairHMM::getMax(double m, double x, double y, unsigned int i, unsigned int j, PairHmmState* state)
{
	if(m >=x && m >=y)
	{
		state->setDiagonalAt(i,j);
		return m;
	}
	else if(x >= y)
	{
		state->setHorizontalAt(i,j);
		return x;
	}
	else
	{
		state->setVerticalAt(i,j);
		return y;
	}

}

void ViterbiPairHMM::runViterbiAlgorithm()
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

			xm = (*M)(i-1,j) + X->getTransitionProbability(M);
			xx = (*X)(i-1,j) + X->getTransitionProbability(X);
			xy = (*Y)(i-1,j) + X->getTransitionProbability(Y);

			X->setValue(i,j,getMax(xm,xx,xy,i,j,X) + emissionX);

			ym = (*M)(i,j-1) + Y->getTransitionProbability(M);
			yx = (*X)(i,j-1) + Y->getTransitionProbability(X);
			yy = (*Y)(i,j-1) + Y->getTransitionProbability(Y);

			Y->setValue(i,j,getMax(ym,yx,yy,i,j,Y) + emissionY);

			mm = (*M)(i-1,j-1) + M->getTransitionProbability(M);
			mx = (*X)(i-1,j-1) + M->getTransitionProbability(X);
			my = (*Y)(i-1,j-1) + M->getTransitionProbability(Y);

			M->setValue(i,j,getMax(mm,mx,my,i,j,M) + emissionM);
		}
	}
}

} /* namespace EBC */
