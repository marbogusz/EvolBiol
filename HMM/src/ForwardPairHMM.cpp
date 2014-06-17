/*
 * ForwardPairHMM.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "Definitions.hpp"
#include "ForwardPairHMM.hpp"
#include "PairwiseHmmInsertState.hpp"
#include "PairwiseHmmDeleteState.hpp"
#include "PairwiseHmmMatchState.hpp"
#include "GTRModel.hpp"
#include "HKY85Model.hpp"
#include "AminoacidSubstitutionModel.hpp"
#include "DpMatrixLoMem.hpp"

namespace EBC
{


ForwardPairHMM::BFGS::BFGS(ForwardPairHMM* enclosing, Definitions::OptimizationType ot) : optimizationType(ot)
{
	parent = enclosing;
	paramsCount = enclosing->optParametersCount;
	this->initParams.set_size(paramsCount);
	this->lowerBounds.set_size(paramsCount);
	this->upperBounds.set_size(paramsCount);
	int i;
	unsigned int ptr =0;

	if(enclosing->estimateSubstitutionParams)
	{
		for (i=0; i<enclosing->substParameters-1; i++)
		{
			initParams(i) = enclosing->optParameters[i];
			//default probs bounds
			lowerBounds(i) = 0.000001;
			//FIXME - provide external bounds????????????
			upperBounds(i) = 5;
		}
		ptr += enclosing->substParameters-1;
	}
	//set time
	initParams(ptr) = enclosing->optParameters[ptr];
	//default probs bounds
	lowerBounds(ptr) = 0.000001;
	//FIXME - provide external bounds????????????
	upperBounds(ptr) = 3;
	ptr++;

	if(enclosing->estimateIndelParams)
	{
		//indel rate FIXME
		initParams(ptr) = enclosing->optParameters[ptr];
		lowerBounds(ptr) = 0.000001;
		upperBounds(ptr) = 0.2;
		//geometric rate
		initParams(ptr+1) = enclosing->optParameters[ptr+1];
		lowerBounds(ptr+1) = 0.000001;
		upperBounds(ptr+1) = 0.999999;

		/*for (i=0; i<enclosing->indelParameters-1; i++)
		{
			initParams(ptr+i) = enclosing->optParameters[ptr+i];
		//default probs bounds
			lowerBounds(ptr+i) = 0.000001;
		//FIXME - provide external bounds????????????
			upperBounds(ptr+i) = 0.999999;
		}*/
	}
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

	switch(optimizationType)
	{
		case Definitions::OptimizationType::BFGS:
		{
			dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
					dlib::objective_delta_stop_strategy(1e-8),
					f_objective,
					derivative(f_objective),
					initParams,
					lowerBounds,
					upperBounds);
			break;
		}
		case Definitions::OptimizationType::BOBYQA:
		{
			dlib::find_min_bobyqa(f_objective, initParams, 10,
					lowerBounds,upperBounds, 0.05, 1e-6, 10000 );
			break;
		}
	}
	DEBUG("BFGS return: " << initParams );
}


ForwardPairHMM::ForwardPairHMM(Sequences* inputSeqs, Definitions::ModelType model ,std::vector<double> indel_params,
		std::vector<double> subst_params, Definitions::OptimizationType ot, bool banding, unsigned int bandPercentage, double evolDistance,
		unsigned int rateCategories, double alpha, bool estimateAlpha) :
		EvolutionaryPairHMM(inputSeqs), userIndelParameters(indel_params), userSubstParameters(subst_params),
		gammaRateCategories(rateCategories), initialAlpha(alpha), estimateAlpha(estimateAlpha)
{
	DEBUG("Creating the model");
	if (model == Definitions::ModelType::GTR)
	{
		substModel = new GTRModel(dict, maths,gammaRateCategories);
	}
	else if (model == Definitions::ModelType::HKY85)
	{
		substModel = new HKY85Model(dict, maths,gammaRateCategories);
	}
	else if (model == Definitions::ModelType::LG)
	{
			substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaLgModel);
	}

	substModel->setAlpha(initialAlpha);

	estimateSubstitutionParams = (subst_params.size() == 0 && substModel->getParamsNumber() > 1);
	estimateIndelParams = (indel_params.size() == 0);
	estimateDivergence = (evolDistance < 0);

	//FIXME
	//Hardcode the band for now
	bandFactor = bandPercentage;
	bandingEnabled = banding;

	//initialize parameter arrays
	initializeModels();


	//TODO - set parameters depending on the values provided
	setParameters();
	if(!estimateDivergence)
	{
		this->mlParameters[this->substParameters-1] = evolDistance;
	}

	getSequencePair();
	getBandWidth();
	calculateModels();
	initializeStates();

	if (estimateDivergence)
	{
		this-> bfgs = new BFGS(this,ot);
		bfgs->optimize();
	}
	else
		this->runForwardAlgorithm();
}

ForwardPairHMM::~ForwardPairHMM()
{
	// TODO Auto-generated destructor stub
	delete bfgs;
	//delete Y;
	//delete X;
	//delete M;
	delete[] optParameters;
	delete[] mlParameters;
	//delete substModel;
	//delete indelModel;
    delete maths;
}

void ForwardPairHMM::initializeStates()
{

	if (M != NULL)
		delete M;
	if (X != NULL)
		delete X;
	if (Y != NULL)
		delete Y;

	//M = new PairwiseHmmMatchState(xSize,ySize);
	//X = new PairwiseHmmInsertState(xSize,ySize);
	//Y = new PairwiseHmmDeleteState(xSize,ySize);

	M = new PairwiseHmmMatchState(new DpMatrixLoMem(xSize,ySize));
	X = new PairwiseHmmInsertState(new DpMatrixLoMem(xSize,ySize));
	Y = new PairwiseHmmDeleteState(new DpMatrixLoMem(xSize,ySize));
}


double ForwardPairHMM::runForwardIteration(const column_vector& bfgsParameters)
{
	unsigned int optPointer = estimateSubstitutionParams ? 0 : substParameters -1;
	cerr << " optimizing : ";
	for(int i=0; i<optParametersCount; i++)
	{
		mlParameters[i+optPointer] = bfgsParameters(i);
		cerr << "\t" << bfgsParameters(i);
	}

	//paste the opt parameters to mlVector
	DEBUGV(mlParameters, totalParameters);

	return this->runForwardAlgorithm() * -1;
}

void ForwardPairHMM::setParameters()
{
	//Models are initialized at this stage
	double* initialSubstData;
	double* initialIndelData;
	double initialDistance = this->generateInitialDistanceParameter();
	unsigned int optPointer = 0;

	if (estimateSubstitutionParams)
	{
		initialSubstData=this->generateInitialSubstitutionParameters();
		for (int i=0; i< this->substParameters-1; i++)
		{
			this->mlParameters[i] = initialSubstData[i];
			this->optParameters[i]= initialSubstData[i];
		}
		delete[] initialSubstData;
		optPointer += substParameters -1;
	}
	else
	{
		for (int i=0; i< this->substParameters-1; i++)
		{
			this->mlParameters[i] = userSubstParameters[i];
		}
	}

	//set time

	this->mlParameters[substParameters-1] = initialDistance;
	this->optParameters[optPointer] = initialDistance;
	optPointer++;

	if(estimateIndelParams)
	{
		initialIndelData = this->generateInitialIndelParameters();
		for (int i=0; i< this->indelParameters-1; i++)
		{
			this->mlParameters[i+substParameters] = initialIndelData[i];
			this->optParameters[i+optPointer]= initialIndelData[i];
		}
		delete[] initialIndelData;
	}
	else
	{
		//set provided vales
		for (int i=0; i< this->indelParameters-1; i++)
		{
			this->mlParameters[i+substParameters] = userIndelParameters[i];
		}
	}

}

double ForwardPairHMM::runForwardAlgorithm()
{

	calculateModels();
	setTransitionProbabilities();

	unsigned int i;
	unsigned int j;
	unsigned int k;
	unsigned int l;

	double sX,sY,sM, sS;

	double xx,xy,xm,yx,yy,ym,mx,my,mm;

	double emissionM;
	double emissionX;
	double emissionY;


	//ReducedPairHmmMatchState* pM = static_cast<ReducedPairHmmMatchState*>(M);
	//ReducedPairHmmInsertState* pX = static_cast<ReducedPairHmmInsertState*>(X);
	//ReducedPairHmmDeleteState* pY = static_cast<ReducedPairHmmDeleteState*>(Y);


	M->initializeData();
	X->initializeData();
	Y->initializeData();


	for (i = 0; i<xSize; i++)
	{
		for (j = 0; j<ySize; j++)
		{
			if(this->withinBand(i,j,this->bandSpan) || !bandingEnabled)
			{
				if(i!=0)
				{
					k = i-1;
					emissionX = log(substModel->getQXi(seq1[i-1].getMatrixIndex()));
					xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
					xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
					xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();
					X->setValueAt(i,j, emissionX + maths->logSum(xm,xx,xy));
				}

				if(j!=0)
				{
					k = j-1;
					emissionY = log(substModel->getQXi(seq2[j-1].getMatrixIndex()));
					ym = M->getValueAt(i,k) + Y->getTransitionProbabilityFromMatch();
					yx = X->getValueAt(i,k) + Y->getTransitionProbabilityFromInsert();
					yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
					Y->setValueAt(i,j, emissionY + maths->logSum(ym,yx,yy));
				}

				if(i!=0 && j!=0 )
				{
					k = i-1;
					l = j-1;
					emissionM = log(substModel->getPXiYi(seq1[i-1].getMatrixIndex(), seq2[j-1].getMatrixIndex()));
					mm = M->getValueAt(k,l) + M->getTransitionProbabilityFromMatch();
					mx = X->getValueAt(k,l) + M->getTransitionProbabilityFromInsert();
					my = Y->getValueAt(k,l) + M->getTransitionProbabilityFromDelete();
					M->setValueAt(i,j, emissionM + maths->logSum(mm,mx,my));
				}
			}
		}
		//pM->outputRow();
	}

	sM = M->getValueAt(xSize-1, ySize-1);
	sX = X->getValueAt(xSize-1, ySize-1);
	sY = Y->getValueAt(xSize-1, ySize-1);
	sS = maths->logSum(sM,sX,sY);

	cerr << "\t" << sX << "\t" << sY << "\t"<< sM << "\t" << sS << endl;

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
	//set the counts and params!
	this->indelParameters = indelModel->getParamsNumber();
	this->substParameters = substModel->getParamsNumber();
	this->totalParameters = indelParameters + substParameters -1;
	this->mlParameters = new double[totalParameters];

	unsigned int substMLParamCount = this->estimateSubstitutionParams == true ? this->substParameters-1 : 0;
	unsigned int indelMLParamCount = this->estimateIndelParams == true ? this->indelParameters-1 : 0;

	optParametersCount = substMLParamCount  + indelMLParamCount + 1;
	//estimate the distance at the minimum
	this->optParameters = new double[optParametersCount];

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
}

} /* namespace EBC */
