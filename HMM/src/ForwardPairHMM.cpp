/*
 * ForwardPairHMM.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "ForwardPairHMM.hpp"


namespace EBC
{


ForwardPairHMM::BFGS::BFGS(ForwardPairHMM* enclosing)
{
	parent = enclosing;
	paramsCount = enclosing->totalParameters;
	this->initParams.set_size(paramsCount);
	this->lowerBounds.set_size(paramsCount);
	this->upperBounds.set_size(paramsCount);

	for (int i=0; i<paramsCount; i++)
	{
		initParams(i) = enclosing->mlParameters[i];
		//default probs bounds
		lowerBounds(i) = 0.000001;
		upperBounds(i) = 1.0;
	}

	//FIXME - hardcoding bounds
	//kappa

	upperBounds(0) = 5;

	DEBUG("++++++++BFGS init with " << paramsCount << " parameters");
}

ForwardPairHMM::BFGS::~BFGS()
{
}

double ForwardPairHMM::BFGS::objectiveFunction(const column_vector& bfgsParameters)
{
	return parent->runForwardIteration(bfgsParameters);
}


const column_vector ForwardPairHMM::BFGS::objectiveFunctionDerivative(const column_vector& bfgsParameters)
{
	column_vector results(this->paramsCount);
	return results;
}


void ForwardPairHMM::BFGS::optimize()
{
	using std::placeholders::_1;
	std::function<double(const column_vector&)> f_objective= std::bind( &ForwardPairHMM::BFGS::objectiveFunction, this, _1 );

/*	dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
			dlib::objective_delta_stop_strategy(1e-9),
			f_objective,
			derivative(f_objective),
			initParams,
			lowerBounds,
			upperBounds);

*/

	dlib::find_min_bobyqa(f_objective, initParams, 10,
			lowerBounds,upperBounds, 0.05, 1e-6, 350 );

	DEBUG("BFGS return: " << initParams );
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

double ForwardPairHMM::runForwardIteration(const column_vector& bfgsParameters)
{
	for(int i=0; i<totalParameters; i++)
	{
		mlParameters[i] = bfgsParameters(i);
	}


	DEBUGV(mlParameters, totalParameters);

	calculateModels();
	initializeStates();
	setTransitionProbabilities();
	//minimize the negative


	//this->substModel->summarize();
	//this->indelModel->summarize();

	return this->runForwardAlgorithm() * -1;
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


	DEBUG ("Forward results:");
	DEBUG (" sX, sY, sM, sS " << sX << "\t" << sY << "\t" << sM << "\t" << sS);

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
