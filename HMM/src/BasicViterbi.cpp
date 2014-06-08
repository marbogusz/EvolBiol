/*
 * BasicViterbi.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 *
 *      An ugly mockup of viterbi!
 *      Change to a common core!!
 */

#include "BasicViterbi.hpp"
#include "GTRModel.hpp"
#include "HKY85Model.hpp"
#include "AminoacidSubstitutionModel.hpp"
#include "AffineGeometricGapModel.hpp"
#include "PairHmmMatchState.hpp"
#include "PairHmmInsertionState.hpp"
#include "PairHmmDeletionState.hpp"
#include "Definitions.hpp"
#include <sstream>

namespace EBC
{

BasicViterbi::BasicViterbi(Sequences* inputSeqs, Definitions::ModelType model,std::vector<double> substParams, double distance, std::vector<double> indelParams, double* estimatedParams)
	: inputSequences(inputSeqs)
{
	//FIXME - deal with alpha!

	M = X = Y = NULL;

	dict = inputSeqs->getDictionary();
	
    maths  = new Maths();
	
    DEBUG("Creating the model");
    if (model == Definitions::ModelType::GTR)
    {
    	substModel = new GTRModel(dict, maths,0);
    }
    else if (model == Definitions::ModelType::HKY85)
    {
    	substModel = new HKY85Model(dict, maths,0);
    }
    else if (model == Definitions::ModelType::LG)
    {
    	substModel = new AminoacidSubstitutionModel(dict, maths,0,Definitions::aaLgModel);
    }

	
    DEBUG("Creating the gap model");
	indelModel = new AffineGeometricGapModel();

	DEBUG("initializeModels");
	initializeModels();

	if (estimatedParams != NULL)
	{
		for (unsigned int idx=0; idx < totalParameters; idx++)
		{
			mlParameters[idx] = estimatedParams[idx];
		}
	}
	else
	{
		//Now copy ml params to the model
		if (totalParameters != substParams.size()+indelParams.size()+1)
		{
			throw ProgramException("Provided parameters count does not match the model\n");
		}
		throw ProgramException("Not implemented yet\n");
	}

	DEBUG("calculateModels");
	calculateModels();
	
    DEBUG("gettingSequences");
	getSequencePair();
	
    DEBUG("initializeStates");
	initializeStates();
	
    DEBUG("setting transition probabilities");

}


void BasicViterbi::initializeModels()
{
	this->indelParameters = indelModel->getParamsNumber();
	this->substParameters = substModel->getParamsNumber();
	this->totalParameters = indelParameters + substParameters -1;
	this->mlParameters = new double[totalParameters];

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
}

void BasicViterbi::initializeStates()
{
	double e,g;
	e = indelModel->getGapExtensionProbability();
	g = indelModel->getGapOpeningProbability();

	if (M != NULL)
		delete M;
	if (X != NULL)
		delete X;
	if (Y != NULL)
		delete Y;

	M = new PairHmmMatchState(xSize,ySize,g,e);
	X = new PairHmmInsertionState(xSize,ySize,g,e);
	Y = new PairHmmDeletionState(xSize,ySize,g,e);

	M->addTransitionProbabilityFrom(M,log(1.0-2.0*g));
	M->addTransitionProbabilityFrom(X,log((1.0-e)*(1.0-2.0*g)));
	M->addTransitionProbabilityFrom(Y,log((1.0-e)*(1.0-2.0*g)));

	X->addTransitionProbabilityFrom(X,log(e+((1-e)*g)));
	Y->addTransitionProbabilityFrom(Y,log(e+((1-e)*g)));

	X->addTransitionProbabilityFrom(Y,log((1.0-e)*g));
	Y->addTransitionProbabilityFrom(X,log((1.0-e)*g));

	X->addTransitionProbabilityFrom(M,log(g));
	Y->addTransitionProbabilityFrom(M,log(g));
}


void BasicViterbi::getSequencePair()
{
	this->seq1 = inputSequences->getSequencesAt(0);
	this->seq2 = inputSequences->getSequencesAt(1);
	this->xSize = seq1.size() +1;
	this->ySize = seq2.size() +1;
}


void BasicViterbi::calculateModels()
{
	indelModel->setParameters(mlParameters+substParameters-1);
	substModel->setParameters(mlParameters);
	//substModel->setParametersInMatrix();
	//substModel->setDiagMeans();
	//substModel->doEigenDecomposition();
	substModel->calculatePt();
}

double BasicViterbi::getMax(double m, double x, double y, unsigned int i, unsigned int j, PairHmmState* state)
{
	if(m >=x && m >=y)
	{
		//state->setDiagonalAt(i,j);

		//cout << i << " " << j << " coming from M" << endl;
		state->setDirection(i,j);
		state->setSrc(i,j,M);
		return m;
	}
	else if(x >= y)
	{

		//state->setVerticalAt(i,j);
		//cout << i << " " << j << " coming from X" << endl;
		state->setDirection(i,j);
		state->setSrc(i,j,X);
		return x;
	}
	else
	{
		//state->setHorizontalAt(i,j);
		//cout << i << " " << j << " coming from Y" << endl;
		state->setDirection(i,j);
		state->setSrc(i,j,Y);
		return y;
	}

}

void BasicViterbi::getResults(stringstream& ss)
{

	double mv, xv, yv;

	//cout << "M" << endl;
	//M->outputValues(0);
	//cout << "X" << endl;
	//X->outputValues(0);
	//cout << "Y" << endl;
	//Y->outputValues(0);

	//DEBUG("Get Results");
	pair<string, string> initialAlignment = std::make_pair("","");
	initialAlignment.first.reserve(xSize > ySize ? xSize : ySize);
	initialAlignment.second.reserve(xSize > ySize ? xSize : ySize);

	string a = inputSequences->getRawSequenceAt(0);
	string b = inputSequences->getRawSequenceAt(1);

	mv =  (*M)(xSize-1,ySize-1) ;
	xv =  (*X)(xSize-1,ySize-1) ;
	yv =  (*Y)(xSize-1,ySize-1) ;


	if(mv >=xv && mv >=yv)
	{
		M->traceback(a,b, &initialAlignment);
	}
	else if(xv >= yv)
	{
		X->traceback(a,b, &initialAlignment);
	}
	else
	{
		Y->traceback(a,b, &initialAlignment);
	}



	//M->traceback(a,b, &initialAlignment);

	reverse(initialAlignment.first.begin(), initialAlignment.first.end());
	reverse(initialAlignment.second.begin(), initialAlignment.second.end());


	ss << ">A" << endl;
	ss << initialAlignment.first << endl;
	ss << ">B" << endl;
	ss << initialAlignment.second << endl;
}

void BasicViterbi::runViterbiAlgorithm()
{
	DEBUG("Run Viterbi");

	unsigned int i;
	unsigned int j;

	double xx,xy,xm,yx,yy,ym,mx,my,mm;

	double emissionM;
	double emissionX;
	double emissionY;

	//while (i != xSize && j != ySize)


	for (i = 0; i<xSize; i++)
	{
		for (j = 0; j<ySize; j++)
		{

			if(i!=0)
			{
				emissionX = log(substModel->getQXi(seq1[i-1].getMatrixIndex()));

				xm = (*M)(i-1,j);
				xm += X->getTransitionProbabilityFrom(M);
				xx = (*X)(i-1,j);
				xx += X->getTransitionProbabilityFrom(X);
				xy = (*Y)(i-1,j);
				xy += X->getTransitionProbabilityFrom(Y);

				//cout << "X state ";
				X->setValue(i,j,getMax(xm,xx,xy,i,j,X) + emissionX);
			}
			if(j!=0)
			{
				emissionY = log(substModel->getQXi(seq2[j-1].getMatrixIndex()));

				ym = (*M)(i,j-1);
				ym += Y->getTransitionProbabilityFrom(M);
				yx = (*X)(i,j-1);
				yx += Y->getTransitionProbabilityFrom(X);
				yy = (*Y)(i,j-1);
				yy += Y->getTransitionProbabilityFrom(Y);
				//cout << "Y state ";
				Y->setValue(i,j,getMax(ym,yx,yy,i,j,Y) + emissionY);
			}

			if(i!=0 && j!=0)
			{
				emissionM = log(substModel->getPXiYi(seq1[i-1].getMatrixIndex(), seq2[j-1].getMatrixIndex()));

				mm = (*M)(i-1,j-1);
				mm += M->getTransitionProbabilityFrom(M);
				mx = (*X)(i-1,j-1);
				mx += M->getTransitionProbabilityFrom(X);
				my = (*Y)(i-1,j-1);
				my += M->getTransitionProbabilityFrom(Y);
				//cout << "M state ";
				M->setValue(i,j,getMax(mm,mx,my,i,j,M) + emissionM);
			}
		}
	}
	DEBUG("Final Viterbi M  " << (*M)(xSize-1,ySize-1) );
	DEBUG("Final Viterbi X  " << (*X)(xSize-1,ySize-1) );
	DEBUG("Final Viterbi Y  " << (*Y)(xSize-1,ySize-1) );
}


BasicViterbi::~BasicViterbi()
{
	// TODO Auto-generated destructor stub
}

} /* namespace EBC */
