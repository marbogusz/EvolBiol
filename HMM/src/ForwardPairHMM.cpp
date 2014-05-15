/*
 * ForwardPairHMM.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "ForwardPairHMM.hpp"
#include "ReducedPairHmmInsertState.hpp"
#include "ReducedPairHmmDeleteState.hpp"
#include "ReducedPairHmmMatchState.hpp"

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

	DEBUG("DLIB optimizer init with " << paramsCount << " parameters");
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
/*
	dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
			dlib::objective_delta_stop_strategy(1e-8),
			f_objective,
			derivative(f_objective),
			initParams,
			lowerBounds,
			upperBounds);

*/
	dlib::find_min_bobyqa(f_objective, initParams, 10,
			lowerBounds,upperBounds, 0.05, 1e-6, 10000 );

	DEBUG("BFGS return: " << initParams );
}


ForwardPairHMM::ForwardPairHMM(Sequences* inputSeqs, bool optimize) :
		EvolutionaryPairHMM(inputSeqs)
{

	bandFactor = 30;

	if (optimize)
	{
		initializeModels();
		getSequencePair();
		getBandWidth();
		calculateModels();
		initializeStates();
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
		getBandWidth();

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

void ForwardPairHMM::initializeStates()
{

	if (M != NULL)
		delete M;
	if (X != NULL)
		delete X;
	if (Y != NULL)
		delete Y;

	M = new ReducedPairHmmMatchState(xSize,ySize);
	X = new ReducedPairHmmInsertState(xSize,ySize);
	Y = new ReducedPairHmmDeleteState(xSize,ySize);
}


double ForwardPairHMM::runForwardIteration(const column_vector& bfgsParameters)
{
	for(int i=0; i<totalParameters; i++)
	{
		mlParameters[i] = bfgsParameters(i);
	}


	DEBUGV(mlParameters, totalParameters);

	return this->runForwardAlgorithm() * -1;
}



double ForwardPairHMM::runForwardAlgorithm()
{

	calculateModels();
	setTransitionProbabilities();

	unsigned int i;
	unsigned int j;

	double sX,sY,sM, sS;

	double xx,xy,xm,yx,yy,ym,mx,my,mm;

	double emissionM;
	double emissionX;
	double emissionY;


	ReducedPairHmmMatchState* pM = static_cast<ReducedPairHmmMatchState*>(M);
	ReducedPairHmmInsertState* pX = static_cast<ReducedPairHmmInsertState*>(X);
	ReducedPairHmmDeleteState* pY = static_cast<ReducedPairHmmDeleteState*>(Y);


	pM->initializeData();
	pX->initializeData();
	pY->initializeData();


	for (i = 0; i<xSize; i++)
	{
		for (j = 0; j<ySize; j++)
		{
			if(this->withinBand(i,j,this->bandSpan))
			{
				if(i!=0)
				{
					emissionX = log(substModel->getQXi(seq1[i-1].getMatrixIndex()));
					xm = pM->valueAtTop(j) + X->getTransitionProbabilityFromMatch();
					xx = pX->valueAtTop(j) + X->getTransitionProbabilityFromInsert();
					xy = pY->valueAtTop(j) + X->getTransitionProbabilityFromDelete();
					pX->setValue(j, emissionX + maths->logSum(xm,xx,xy));
				}

				if(j!=0)
				{
					emissionY = log(substModel->getQXi(seq2[j-1].getMatrixIndex()));
					ym = pM->valueAtLeft(j) + Y->getTransitionProbabilityFromMatch();
					yx = pX->valueAtLeft(j) + Y->getTransitionProbabilityFromInsert();
					yy = pY->valueAtLeft(j) + Y->getTransitionProbabilityFromDelete();
					pY->setValue(j, emissionY + maths->logSum(ym,yx,yy));
				}

				if(i!=0 && j!=0 )
				{
					emissionM = log(substModel->getPXiYi(seq1[i-1].getMatrixIndex(), seq2[j-1].getMatrixIndex()));
					mm = pM->valueAtDiagonal(j) + M->getTransitionProbabilityFromMatch();
					mx = pX->valueAtDiagonal(j) + M->getTransitionProbabilityFromInsert();
					my = pY->valueAtDiagonal(j) + M->getTransitionProbabilityFromDelete();
					pM->setValue(j, emissionM + maths->logSum(mm,mx,my));
				}
			}
		}

		//pM->outputRow();
		pX->nextRow();
		pY->nextRow();
		pM->nextRow();
	}



	pX->nextRow();
	pY->nextRow();
	pM->nextRow();

	sM = pM->valueAtColumn(ySize-1);
	sX = pX->valueAtColumn(ySize-1);
	sY = pY->valueAtColumn(ySize-1);
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

void ForwardPairHMM::initializeModels()
{
	generateInitialParameters();
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

} /* namespace EBC */
