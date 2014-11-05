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


BackwardPairHMM::BackwardPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, SubstitutionModelBase* smdl,
		IndelModel* imdl,  Definitions::DpMatrixType mt ,Band* bandObj) :
		EvolutionaryPairHMM(s1,s2, smdl, imdl, mt, bandObj)
{
}

BackwardPairHMM::~BackwardPairHMM()
{
}

void BackwardPairHMM::calculatePosteriors(ForwardPairHMM* fwd)
{
	int i,j;
	double xval, yval, mval;
	double fwdM;

	fwdM = fwd->M->getValueAt(xSize-1,ySize-1);

	for (i = 0; i<xSize-1; i++)
	{
		for (j = 0; j<ySize-1; j++)
		{
			xval = X->getValueAt(i,j) + fwd->X->getValueAt(i,j) - fwdM;
			yval = Y->getValueAt(i,j) + fwd->Y->getValueAt(i,j) - fwdM;
			mval = M->getValueAt(i,j) + fwd->M->getValueAt(i,j) - fwdM;

			X->setValueAt(i,j,xval);
			Y->setValueAt(i,j,yval);
			M->setValueAt(i,j,mval);
		}
	}

}

double BackwardPairHMM::runAlgorithm()
{
	calculateModels();
	setTransitionProbabilities();

	int i;
	int j;

	double sX,sY,sM, sS;

	double bm, bx, by;

	double bmp, bxp,byp;

	M->initializeData(true);
	X->initializeData(true);
	Y->initializeData(true);

	for (j = ySize-1, i=xSize-1; j >= 0; j--)
	{
		bxp = (i==xSize-1) ? -100000 : X->getValueAt(i+1,j) + ptmatrix->getLogEquilibriumFreq(seq1[i].getMatrixIndex());
		byp = (j==ySize-1) ? -100000 : Y->getValueAt(i,j+1) + ptmatrix->getLogEquilibriumFreq(seq2[j].getMatrixIndex());
		bmp = (i==xSize-1 ||j==ySize-1) ? -100000 : M->getValueAt(i+1,j+1) + ptmatrix->getLogPairTransition(seq1[i].getMatrixIndex(), seq2[j].getMatrixIndex());

		bx = maths->logSum(M->getTransitionProbabilityFromInsert() +  bmp,
				X->getTransitionProbabilityFromInsert() + bxp,
				Y->getTransitionProbabilityFromInsert() + byp);

		by = maths->logSum(M->getTransitionProbabilityFromDelete() + bmp,
				X->getTransitionProbabilityFromDelete() + bxp,
				Y->getTransitionProbabilityFromDelete() + byp);

		bm = maths->logSum(M->getTransitionProbabilityFromMatch() + bmp,
				X->getTransitionProbabilityFromMatch() + bxp,
				Y->getTransitionProbabilityFromMatch() + byp);

		if (i==xSize-1  && j==ySize-1)
		{
			X->setValueAt(i, j, 0);
			Y->setValueAt(i, j, 0);
			M->setValueAt(i, j, 0);
		}
		else
		{
			X->setValueAt(i, j, bx);
			Y->setValueAt(i, j, by);
			M->setValueAt(i, j, bm);
		}
	}
	for (i = xSize-1, j=ySize-1; i >= 0; i--)
	{
		bxp = (i==xSize-1) ? -100000 : X->getValueAt(i+1,j) + ptmatrix->getLogEquilibriumFreq(seq1[i].getMatrixIndex());
		byp = (j==ySize-1) ? -100000 : Y->getValueAt(i,j+1) + ptmatrix->getLogEquilibriumFreq(seq2[j].getMatrixIndex());
		bmp = (i==xSize-1 ||j==ySize-1) ? -100000 : M->getValueAt(i+1,j+1) + ptmatrix->getLogPairTransition(seq1[i].getMatrixIndex(), seq2[j].getMatrixIndex());

		bx = maths->logSum(M->getTransitionProbabilityFromInsert() + bmp,
				X->getTransitionProbabilityFromInsert() + bxp,
				Y->getTransitionProbabilityFromInsert() + byp);

		by = maths->logSum(M->getTransitionProbabilityFromDelete() + bmp,
				X->getTransitionProbabilityFromDelete() + bxp,
				Y->getTransitionProbabilityFromDelete() + byp);

		bm = maths->logSum(M->getTransitionProbabilityFromMatch() + bmp,
				X->getTransitionProbabilityFromMatch() + bxp,
				Y->getTransitionProbabilityFromMatch() + byp);

		if (i==xSize-1  && j==ySize-1)
		{
			X->setValueAt(i, j, 0);
			Y->setValueAt(i, j, 0);
			M->setValueAt(i, j, 0);
		}
		else
		{
			X->setValueAt(i, j, bx);
			Y->setValueAt(i, j, by);
			M->setValueAt(i, j, bm);
		}
	}


	for (i = xSize-2; i >= 0; i--)
	{
		for (j = ySize-2; j >= 0; j--)
		{

			bxp = X->getValueAt(i+1,j) + ptmatrix->getLogEquilibriumFreq(seq1[i].getMatrixIndex());
			byp = Y->getValueAt(i,j+1) + ptmatrix->getLogEquilibriumFreq(seq2[j].getMatrixIndex());
			bmp = M->getValueAt(i+1,j+1) + ptmatrix->getLogPairTransition(seq1[i].getMatrixIndex(), seq2[j].getMatrixIndex());

			bx = maths->logSum(M->getTransitionProbabilityFromInsert() + bmp,
					X->getTransitionProbabilityFromInsert() + bxp,
					Y->getTransitionProbabilityFromInsert() + byp);

			by = maths->logSum(M->getTransitionProbabilityFromDelete() + bmp,
								X->getTransitionProbabilityFromDelete() + bxp,
								Y->getTransitionProbabilityFromDelete() + byp);

			bm = maths->logSum(M->getTransitionProbabilityFromMatch() + bmp,
											X->getTransitionProbabilityFromMatch() + bxp,
											Y->getTransitionProbabilityFromMatch() + byp);

				X->setValueAt(i, j, bx);
				Y->setValueAt(i, j, by);
				M->setValueAt(i, j, bm);

		}
	}

	sM = M->getValueAt(0, 0);
	sX = X->getValueAt(0, 0);
	sY = Y->getValueAt(0, 0);
	sS = maths->logSum(sM,sX,sY);

	DEBUG ("Backward results:");
	DEBUG (" sX, sY, sM, sS " << sX << "\t" << sY << "\t" << sM << "\t" << sS);

	return sS* -1.0;
}


} /* namespace EBC */
