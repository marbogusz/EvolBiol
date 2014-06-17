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
#include "PairwiseHmmMatchState.hpp"
#include "PairwiseHmmInsertState.hpp"
#include "PairwiseHmmDeleteState.hpp"
#include "Definitions.hpp"
#include <sstream>

namespace EBC
{

BasicViterbi::BasicViterbi(Sequences* inputSeqs, Definitions::ModelType model,std::vector<double> substParams, double distance,
		std::vector<double> indelParams, unsigned int rateCategories, double alpha, double* estimatedParams)
	: inputSequences(inputSeqs)
{
	//FIXME - deal with alpha!

	M = X = Y = NULL;

	dict = inputSeqs->getDictionary();
	
    maths  = new Maths();
	
    DEBUG("Creating the model");
    if (model == Definitions::ModelType::GTR)
    {
    	substModel = new GTRModel(dict, maths,rateCategories);
    }
    else if (model == Definitions::ModelType::HKY85)
    {
    	substModel = new HKY85Model(dict, maths,rateCategories);
    }
    else if (model == Definitions::ModelType::LG)
    {
    	substModel = new AminoacidSubstitutionModel(dict, maths,rateCategories,Definitions::aaLgModel);
    }

    substModel->setAlpha(alpha);

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
			throw HmmException("Provided parameters count does not match the model\n");
		}
		throw HmmException("Not implemented yet\n");
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

	M = new PairwiseHmmMatchState(xSize,ySize);
	X = new PairwiseHmmInsertState(xSize,ySize);
	Y = new PairwiseHmmDeleteState(xSize,ySize);

	M->setTransitionProbabilityFromMatch(log(1-2*g));
	M->setTransitionProbabilityFromInsert(log((1-e)*(1-2*g)));
	M->setTransitionProbabilityFromDelete(log((1-e)*(1-2*g)));

	X->setTransitionProbabilityFromInsert(log(e+((1-e)*g)));
	Y->setTransitionProbabilityFromDelete(log(e+((1-e)*g)));

	X->setTransitionProbabilityFromDelete(log((1-e)*g));
	Y->setTransitionProbabilityFromInsert(log((1-e)*g));

	X->setTransitionProbabilityFromMatch(log(g));
	Y->setTransitionProbabilityFromMatch(log(g));
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

double BasicViterbi::getMax(double m, double x, double y, unsigned int i, unsigned int j, PairwiseHmmStateBase* state)
{
	if(m >=x && m >=y)
	{
		//state->setDiagonalAt(i,j);

		//cout << i << " " << j << " coming from M" << endl;
		state->setDirection(i,j);
		state->setSourceMatrixPtr(i,j,M);
		return m;
	}
	else if(x >= y)
	{

		//state->setVerticalAt(i,j);
		//cout << i << " " << j << " coming from X" << endl;
		state->setDirection(i,j);
		state->setSourceMatrixPtr(i,j,X);
		return x;
	}
	else
	{
		//state->setHorizontalAt(i,j);
		//cout << i << " " << j << " coming from Y" << endl;
		state->setDirection(i,j);
		state->setSourceMatrixPtr(i,j,Y);
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

	mv =  M->getValueAt(xSize-1,ySize-1) ;
	xv =  X->getValueAt(xSize-1,ySize-1) ;
	yv =  Y->getValueAt(xSize-1,ySize-1) ;


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

	unsigned int i,j,k,l;

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

				k = i-1;
				emissionX = log(substModel->getQXi(seq1[i-1].getMatrixIndex()));
				xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
				xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
				xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();

				X->setValueAt(i,j,getMax(xm,xx,xy,i,j,X) + emissionX);
			}
			if(j!=0)
			{
				k = j-1;
				emissionY = log(substModel->getQXi(seq2[j-1].getMatrixIndex()));
				ym = M->getValueAt(i,k) + Y->getTransitionProbabilityFromMatch();
				yx = X->getValueAt(i,k) + Y->getTransitionProbabilityFromInsert();
				yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
				Y->setValueAt(i,j,getMax(ym,yx,yy,i,j,Y) + emissionY);
			}

			if(i!=0 && j!=0)
			{
				k = i-1;
				l = j-1;
				emissionM = log(substModel->getPXiYi(seq1[i-1].getMatrixIndex(), seq2[j-1].getMatrixIndex()));
				mm = M->getValueAt(k,l) + M->getTransitionProbabilityFromMatch();
				mx = X->getValueAt(k,l) + M->getTransitionProbabilityFromInsert();
				my = Y->getValueAt(k,l) + M->getTransitionProbabilityFromDelete();
				M->setValueAt(i,j,getMax(mm,mx,my,i,j,M) + emissionM);
			}
		}
	}
	DEBUG("Final Viterbi M  " << M->getValueAt(xSize-1,ySize-1) );
	DEBUG("Final Viterbi X  " << X->getValueAt(xSize-1,ySize-1) );
	DEBUG("Final Viterbi Y  " << Y->getValueAt(xSize-1,ySize-1) );
}


BasicViterbi::~BasicViterbi()
{
	// TODO Auto-generated destructor stub
}

} /* namespace EBC */
