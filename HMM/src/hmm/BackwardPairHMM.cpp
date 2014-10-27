/*
 * BackwardPairHMM.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "core/Definitions.hpp"
#include "hmm/BackwardPairHMM.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "hmm/DpMatrixFull.hpp"


namespace EBC
{


BackwardPairHMM::BackwardPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, bool banding,
		SubstitutionModelBase* smdl, IndelModel* imdl, unsigned int bandPercentage, Definitions::DpMatrixType mt) :
		EvolutionaryPairHMM(s1,s2, banding, smdl, imdl, bandPercentage,mt)
{
}

BackwardPairHMM::~BackwardPairHMM()
{
}

void BackwardPairHMM::calculatePosteriors(ForwardPairHMM* fwd)
{
	int i,j;
	double xval, yval, mval;
	double fwdX, fwdY, fwdM;

	DEBUG("Backward Match");
	dynamic_cast<DpMatrixFull*>(M->getDpMatrix())->outputValues(0);

	DEBUG("Backward Insert");
	dynamic_cast<DpMatrixFull*>(X->getDpMatrix())->outputValues(0);

	DEBUG("Backward Delete");
	dynamic_cast<DpMatrixFull*>(Y->getDpMatrix())->outputValues(0);

	DEBUG("Forward Match");
	dynamic_cast<DpMatrixFull*>(fwd->M->getDpMatrix())->outputValues(0);

	DEBUG("Forward Insert");
	dynamic_cast<DpMatrixFull*>(fwd->X->getDpMatrix())->outputValues(0);

	DEBUG("Forward Delete");
	dynamic_cast<DpMatrixFull*>(fwd->Y->getDpMatrix())->outputValues(0);


	fwdX = fwd->X->getValueAt(xSize,ySize);
	fwdY = fwd->Y->getValueAt(xSize,ySize);
	fwdM = fwd->M->getValueAt(xSize,ySize);

	for (i = 1; i<xSize; i++)
	{
		for (j = 1; j<ySize; j++)
		{
			xval = X->getValueAt(i,j) + fwd->X->getValueAt(i,j) - fwdX;
			yval = Y->getValueAt(i,j) + fwd->Y->getValueAt(i,j) - fwdY;
			mval = M->getValueAt(i,j) + fwd->M->getValueAt(i,j) - fwdM;

			X->setValueAt(i,j,xval);
			Y->setValueAt(i,j,xval);
			M->setValueAt(i,j,xval);
		}
	}

	DEBUG("Posterior Match");
	dynamic_cast<DpMatrixFull*>(M->getDpMatrix())->outputValues(0);

	DEBUG("Posterior Insert");
	dynamic_cast<DpMatrixFull*>(X->getDpMatrix())->outputValues(0);

	DEBUG("Posterior Delete");
	dynamic_cast<DpMatrixFull*>(Y->getDpMatrix())->outputValues(0);

}

double BackwardPairHMM::runAlgorithm()
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

	M->initializeData();
	X->initializeData();
	Y->initializeData();

	unsigned int maxXsize = xSize-1;
	unsigned int maxYsize = ySize-1;

	//Reverse sequences = effectively forward algorithm
	for (i = 0; i < xSize; i++)
	{
		for (j = 0; j < xSize; j++)
		{
				if(i!=0)
				{
					k = i-1;
					emissionX = log(ptmatrix->getEquilibriumFreq(seq1[maxXsize - i].getMatrixIndex()));
					xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
					xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
					xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();
					X->setValueAt(i,j, emissionX + maths->logSum(xm,xx,xy));
				}

				if(j!=0)
				{
					k = j-1;
					emissionY = log(ptmatrix->getEquilibriumFreq(seq2[j-1].getMatrixIndex()));
					ym = M->getValueAt(i,k) + Y->getTransitionProbabilityFromMatch();
					yx = X->getValueAt(i,k) + Y->getTransitionProbabilityFromInsert();
					yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
					Y->setValueAt(i,j, emissionY + maths->logSum(ym,yx,yy));
				}

				if(i!=0 && j!=0 )
				{
					k = i-1;
					l = j-1;
					emissionM = log(ptmatrix->getPairTransition(seq1[i-1].getMatrixIndex(), seq2[j-1].getMatrixIndex()));
					mm = M->getValueAt(k,l) + M->getTransitionProbabilityFromMatch();
					mx = X->getValueAt(k,l) + M->getTransitionProbabilityFromInsert();
					my = Y->getValueAt(k,l) + M->getTransitionProbabilityFromDelete();
					M->setValueAt(i,j, emissionM + maths->logSum(mm,mx,my));
				}
		}
	}

	sM = M->getValueAt(0, 0);
	sX = X->getValueAt(0, 0);
	sY = Y->getValueAt(0, 0);
	sS = maths->logSum(sM,sX,sY);

	//cerr << "\t" << sX << "\t" << sY << "\t"<< sM << "\t" << sS << endl;

	DEBUG ("Backward results:");
	DEBUG (" sX, sY, sM, sS " << sX << "\t" << sY << "\t" << sM << "\t" << sS);

	return sS* -1.0;
}


} /* namespace EBC */
