/*
 * ForwardPairHMM.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "ForwardPairHMM.hpp"

namespace EBC
{

lbfgsfloatval_t ForwardPairHMM::BFGS::calculateGradients(const lbfgsfloatval_t* x,
		lbfgsfloatval_t* g, const int n, const lbfgsfloatval_t step)
{

	for (int i=0; i<5; i++)
	{
		cerr << exp(x[i]) << "\t";
	}
	cerr << endl;

	double localLikelihood = parent->runForwardIteration((const double*)x)*-1.0;
	double differentialLikelihood;

	for (int i=0; i<n; i++)
	{
		copyWithSmallDiff(i,x);
		differentialLikelihood = parent->runForwardIteration((const double*) tempParams)*-1.0;
		g[i] = (differentialLikelihood - localLikelihood) * smallDiffInverse;
	}

	return localLikelihood;
}

void ForwardPairHMM::BFGS::copyWithSmallDiff(int pos, const lbfgsfloatval_t* x)
{
	for (int i=0; i<paramsCount; i++)
	{
		tempParams[i] = x[i];
		if (i==pos)
		{
			tempParams[i] += smallDiff;
		}
	}
}

int ForwardPairHMM::BFGS::reportProgress(const lbfgsfloatval_t* x,
		const lbfgsfloatval_t* g, const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step, int n, int k, int ls)
{
	DEBUG("BFGS value:" << fx << " iteration: " << k);
	//for (int i=0; i<5; i++)
	//{
	//	cerr << exp(x[i]) << "\t";
	//}
	//cerr << endl;
	return 0;
}

double ForwardPairHMM::BFGS::logistic(double x)
{
	return 1.0/(1+exp(x*-1));
}

ForwardPairHMM::BFGS::BFGS(ForwardPairHMM* enclosing)
{
	parent = enclosing;
	//get a pointer to params
	this->currentParameters = enclosing->mlParameters;
	//get a local parmas count
	paramsCount = enclosing->totalParameters;
	DEBUG("++++++++libBFGS init with " << paramsCount << " parameters");
	smallDiff = 1e-6;
	smallDiffExp = exp(smallDiff);
	smallDiffInverse = 1000000;
	smallDiffExpInverse = 1/smallDiffExp;

}

ForwardPairHMM::BFGS::~BFGS()
{
}

void ForwardPairHMM::BFGS::optimize()
{
	lbfgs_parameter_t param;
	lbfgs_parameter_init(&param);
	param.max_iterations  = 100;
	//param.wolfe = 0.11;
	//param.gtol = 0.11;
	//param.min_step = 0.001;
	//param.max_step = 1;
	int ret = lbfgs(paramsCount, currentParameters, &likelihood, _evaluate, _progress, this, &param);
	DEBUG("BFGS return: " << ret << " likelihood " << likelihood);
	DEBUGV(currentParameters,paramsCount);
}


ForwardPairHMM::ForwardPairHMM(Sequences* inputSeqs) :
		EvolutionaryPairHMM(inputSeqs)
{
	initializeModels();
	getSequencePair();
	this-> bfgs = new BFGS(this);
	bfgs->optimize();
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
double ForwardPairHMM::runForwardIteration(const double * bfgsParameters)
{
	this->mlParameters = (double*) bfgsParameters;
	calculateModels();
	initializeStates();
	setTransitionProbabilities();
	return this->runForwardAlgorithm();
}


double ForwardPairHMM::runForwardAlgorithm()
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
			emissionM = log(substModel->getPXiYi(seq1[i-1].getMatrixIndex(), seq2[j-1].getMatrixIndex()));
			emissionX = log(substModel->getQXi(seq1[i-1].getMatrixIndex()));
			emissionY = log(substModel->getQXi(seq2[j-1].getMatrixIndex()));

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


			X->setValue(i,j, emissionX + maths->logSum(xm,xx,xy));
			Y->setValue(i,j, emissionY + maths->logSum(ym,yx,yy));
			M->setValue(i,j, emissionM + maths->logSum(mm,mx,my));

		}
	}
	return (*M)(xSize-1,ySize-1);
}

} /* namespace EBC */
