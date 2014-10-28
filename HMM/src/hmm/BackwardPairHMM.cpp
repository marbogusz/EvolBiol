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


	fwdX = fwd->X->getValueAt(xSize-1,ySize-1);
	fwdY = fwd->Y->getValueAt(xSize-1,ySize-1);
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

	int i;
	int j;

	double sX,sY,sM, sS;

	double bm, bx, by;

	double bmp, bxp,byp;

	double emissionM;
	double emissionX;
	double emissionY;

	M->initializeData(true);
	X->initializeData(true);
	Y->initializeData(true);

	//FIXME - 2 potential errors last corner set to 0 or i,j >=0 conditions
	for (i = xSize-2; i >= 0; i--)
	{
		for (j = ySize-2; j >= 0; j--)
		{

			emissionX = log(ptmatrix->getEquilibriumFreq(seq1[i].getMatrixIndex()));
			emissionY = log(ptmatrix->getEquilibriumFreq(seq2[j].getMatrixIndex()));
			emissionM = log(ptmatrix->getPairTransition(seq1[i].getMatrixIndex(), seq2[j].getMatrixIndex()));

			bxp = X->getValueAt(i+1,j);
			byp = Y->getValueAt(i,j+1);
			bmp = M->getValueAt(i+1,j+1);

			bx = maths->logSum(M->getTransitionProbabilityFromInsert() + emissionM + bmp,
					X->getTransitionProbabilityFromInsert() + emissionX + bxp,
					Y->getTransitionProbabilityFromInsert() + emissionY + byp);

			//TODO bx and by should be identical - same for forward.
			//FIXME - speedup by using only 1 indel calculation - see above
			by = maths->logSum(M->getTransitionProbabilityFromDelete() + emissionM + bmp,
								X->getTransitionProbabilityFromDelete() + emissionX + bxp,
								Y->getTransitionProbabilityFromDelete() + emissionY + byp);

			bm = maths->logSum(M->getTransitionProbabilityFromMatch() + emissionM + bmp,
											X->getTransitionProbabilityFromMatch() + emissionX + bxp,
											Y->getTransitionProbabilityFromMatch() + emissionY + byp);

			X->setValueAt(i, j, bx);
			Y->setValueAt(i, j, by);
			M->setValueAt(i, j, bm);
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
