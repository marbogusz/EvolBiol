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
	int i;
	double tmp;

	for (i=0; i<paramsCount; i++)
	{
		//kappa hack
		tmp = Maths::logistic(x[i]);
		tempParams[i] = i==0 ? 3*tmp : tmp;
	}

	double localLikelihood = parent->runForwardIteration((const double*)tempParams)*-1.0;
	double differentialLikelihood;

	for (i=0; i<n; i++)
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
		tempParams[i] = i==0 ? 3*Maths::logistic(x[i]) : Maths::logistic(x[i]);
		if (i==pos)
		{
			tempParams[i] += smallDiff;
		}

		/*tempParams[i] = x[i];
		if (i==pos)
		{
			tempParams[i] += smallDiff;
		}
		tempParams[i] = Maths::logistic(tempParams[i]);
		*/
	}
}

int ForwardPairHMM::BFGS::reportProgress(const lbfgsfloatval_t* x,
		const lbfgsfloatval_t* g, const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step, int n, int k, int ls)
{
	DEBUG("BFGS value:" << fx << " iteration: " << k);
	for (int i=0; i<n; i++)
	{
		cerr << Maths::logistic(x[i]) << "\t";
	}
	cerr << endl;
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


ForwardPairHMM::ForwardPairHMM(Sequences* inputSeqs, bool optimize) :
		EvolutionaryPairHMM(inputSeqs)
{
	if (optimize)
	{
		initializeModels();
		getSequencePair();
		this-> bfgs = new BFGS(this);
		bfgs->optimize();
	}
	else
	{
		this->indelParameters = indelModel->getParamsNumber();
		this->substParameters = substModel->getParamsNumber();
		this->totalParameters = indelParameters + substParameters -1;
		this->mlParameters = new double[totalParameters];

		testFreqs[0] = 0.4;
		testFreqs[1] = 0.2;
		testFreqs[2] = 0.1;
		testFreqs[3] = 0.3;

		mlParameters[0] = 2;
		mlParameters[1] = 0.1;
		mlParameters[2] = 0.05;
		mlParameters[3] = 0.25;

		substModel->setObservedFrequencies(testFreqs);

		getSequencePair();

		calculateModels();
		initializeStates();
		setTransitionProbabilities();
		runForwardAlgorithm();
	}

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

/*
void ForwardPairHMM::initializeModels()
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
	mlParameters[5] = 0.0905152;	//time
	mlParameters[6] = 0.111341;		//lambda
	mlParameters[7] = 0.521122;		//extension prob



		testFreqs[0] = 0.4;
		testFreqs[1] = 0.2;
		testFreqs[2] = 0.1;
		testFreqs[3] = 0.3;


		mlParameters[0] = 2;
		mlParameters[1] = 0.1;
		mlParameters[2] = 0.05;
		mlParameters[3] = 0.5;

		//start time is the first parameter

		//substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
		substModel->setObservedFrequencies(testFreqs);
}
*/
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

	//calculateModels();
	//initializeStates();
	//setTransitionProbabilities();


	unsigned int i;
	unsigned int j;

	double sX,sY,sM, sS;

	double xx,xy,xm,yx,yy,ym,mx,my,mm;

	double emissionM;
	double emissionX;
	double emissionY;

	for (i = 0; i<xSize; i++)
	{
		for (j = 0; j<ySize; j++)
		{
			if(i!=0)
			{
				emissionX = log(substModel->getQXi(seq1[i-1].getMatrixIndex()));
				xm = (*M)(i-1,j) + X->getTransitionProbabilityFrom(M);
				xx = (*X)(i-1,j) + X->getTransitionProbabilityFrom(X);
				xy = (*Y)(i-1,j) + X->getTransitionProbabilityFrom(Y);
				X->setValue(i,j, emissionX + maths->logSum(xm,xx,xy));
			}

			if(j!=0)
			{
				emissionY = log(substModel->getQXi(seq2[j-1].getMatrixIndex()));
				ym = (*M)(i,j-1) + Y->getTransitionProbabilityFrom(M);
				yx = (*X)(i,j-1) + Y->getTransitionProbabilityFrom(X);
				yy = (*Y)(i,j-1) + Y->getTransitionProbabilityFrom(Y);
				Y->setValue(i,j, emissionY + maths->logSum(ym,yx,yy));
			}

			if(i!=0 && j!=0 )
			{
				emissionM = log(substModel->getPXiYi(seq1[i-1].getMatrixIndex(), seq2[j-1].getMatrixIndex()));
				mm = (*M)(i-1,j-1) + M->getTransitionProbabilityFrom(M);
				mx = (*X)(i-1,j-1) + M->getTransitionProbabilityFrom(X);
				my = (*Y)(i-1,j-1) + M->getTransitionProbabilityFrom(Y);
				M->setValue(i,j, emissionM + maths->logSum(mm,mx,my));
			}
		}
	}
	sM = (*M)(xSize-1,ySize-1);
	sX = (*X)(xSize-1,ySize-1);
	sY = (*Y)(xSize-1,ySize-1);
	sS = maths->logSum(sM,sX,sY);


	//DEBUG ("Forward results:");
	//DEBUG (" sX, sY, sM, sS " << sX << "\t" << sY << "\t" << sM << "\t" << sS);

	/*cout << "M" << endl;
	M->outputValues(0);
	cout << "X" << endl;
	X->outputValues(0);
	cout << "Y" << endl;
	Y->outputValues(0);
	 */
	return sS;
}

} /* namespace EBC */
