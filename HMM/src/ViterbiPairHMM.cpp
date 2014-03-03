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


/*
	--- PairHMM model details ---
	Time: t = 0.0905152
	Model: HKY
	IndelRate:   0.111341
	FractionIns: 0.999001
	ProbExtend: 0.521122
	StateFreq:   0.321661 0.302651 0.107554 0.268134
	PT
	0.307973	0.000988417	0.00391426	0.000775816
	0.00111648	0.274434	0.000124826	0.023757
	0.0350102	0.000988417	0.095185	0.000775816
	0.00111648	0.0302673	0.000124826	0.240077
	IndelProb:   0.0100174
	Transition probabilities for HMM: [0] To Match; [1] To Insert; [2] To Delete
	FromMatch   0.979965 0.0100174 0.0100174
	FromIns     0.469284 0.525919 0.00479709
	FromDel     0.469284 0.00479709 0.525919
	Final likelihood estimate: -1760.05
	Forward accurate esimate:  -1760.05
	And the viterbi: 2  1020
*/


	mlParameters[0] = 0.912374;
	mlParameters[1] = 0.051834;
	mlParameters[2] = 0.000010;
	mlParameters[3] = 0.025448;
	mlParameters[4] = 0.000010;
	mlParameters[5] = 0.0905152;		//time
	mlParameters[6] = 0.111341;		//lambda
	mlParameters[7] = 0.521122;		//extension prob

	testFreqs[0] = 0.321661;
	testFreqs[1] = 0.302651;
	testFreqs[2] = 0.107554;
	testFreqs[3] = 0.268134;

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
			emissionM = log(substModel->getPXiYi(seq1[i-1].getMatrixIndex(), seq2[j-1].getMatrixIndex()));
			emissionX = log(substModel->getQXi(seq1[i-1].getMatrixIndex()));
			emissionY = log(substModel->getQXi(seq2[j-1].getMatrixIndex()));

			xm = (*M)(i-1,j);
			xm += X->getTransitionProbabilityFrom(M);
			xx = (*X)(i-1,j);
			xx += X->getTransitionProbabilityFrom(X);
			xy = (*Y)(i-1,j);
			xy += X->getTransitionProbabilityFrom(Y);

			X->setValue(i,j,getMax(xm,xx,xy,i,j,X) + emissionX);

			ym = (*M)(i,j-1);
			ym += Y->getTransitionProbabilityFrom(M);
			yx = (*X)(i,j-1);
			yx += Y->getTransitionProbabilityFrom(X);
			yy = (*Y)(i,j-1);
			yy += Y->getTransitionProbabilityFrom(Y);

			Y->setValue(i,j,getMax(ym,yx,yy,i,j,Y) + emissionY);

			mm = (*M)(i-1,j-1);
			mm += M->getTransitionProbabilityFrom(M);
			mx = (*X)(i-1,j-1);
			mx += M->getTransitionProbabilityFrom(X);
			my = (*Y)(i-1,j-1);
			my += M->getTransitionProbabilityFrom(Y);

			M->setValue(i,j,getMax(mm,mx,my,i,j,M) + emissionM);
		}
	}
	DEBUG("Final Viterbi M  " << (*M)(xSize-1,ySize-1) );
	DEBUG("Final Viterbi X  " << (*X)(xSize-1,ySize-1) );
	DEBUG("Final Viterbi Y  " << (*Y)(xSize-1,ySize-1) );
}

} /* namespace EBC */
