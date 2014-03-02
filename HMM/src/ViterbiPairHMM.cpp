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
	DEBUG("initializeModels");
	initializeModels();
	//initial model calculation
	DEBUG("calculateModels");
	calculateModels();
	DEBUG("gettingSequences");
	getSequencePair();
	DEBUG("initializeStates");
	initializeStates();
	DEBUG("setting transition probabilities");
	//no need to set them again, initialized below
	//setTransitionProbabilities();
}

ViterbiPairHMM::~ViterbiPairHMM()
{
	// TODO Auto-generated destructor stub
}

void ViterbiPairHMM::initializeModels()
{
	//generateInitialParameters();

	 //time is a parameter with both indel and subst, we use 1 common time
	this->indelParameters = indelModel->getParamsNumber();
	this->substParameters = substModel->getParamsNumber();
	this->totalParameters = indelParameters + substParameters -1;
	this->mlParameters = new double[totalParameters];

	mlParameters[0] = 0.912374;
	mlParameters[1] = 0.051834;
	mlParameters[2] = 0.000010;
	mlParameters[3] = 0.025448;
	mlParameters[4] = 0.000010;
	mlParameters[5] = 0.1;		//time
	mlParameters[6] = 0.05;		//lambda
	mlParameters[7] = 0.4;		//extension prob

	testFreqs[0] = 0.26089;
	testFreqs[1] = 0.32737;
	testFreqs[2] = 0.30782;
	testFreqs[3] = 0.10391;

	//start time is the first parameter

	substModel->setObservedFrequencies(testFreqs);
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

void ViterbiPairHMM::getResults()
{
	//M->outputTrace();

	DEBUG("Get Results");
	pair<string, string> initialAlignment = std::make_pair("","");
	initialAlignment.first.reserve(xSize > ySize ? xSize : ySize);
	initialAlignment.second.reserve(xSize > ySize ? xSize : ySize);

	string a = inputSequences->getRawSequenceAt(0);
	string b = inputSequences->getRawSequenceAt(1);


	M->traceback(a,b, &initialAlignment);

	reverse(initialAlignment.first.begin(), initialAlignment.first.end());
	reverse(initialAlignment.second.begin(), initialAlignment.second.end());

	cout << "Resulting alignment : " << endl;
	cout << initialAlignment.first << endl;
	cout << initialAlignment.second << endl;
}

void ViterbiPairHMM::runViterbiAlgorithm()
{
	DEBUG("Run Viterbi");

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
