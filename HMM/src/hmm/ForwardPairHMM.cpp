/*
 * ForwardPairHMM.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "core/Definitions.hpp"
#include "hmm/ForwardPairHMM.hpp"

namespace EBC
{


ForwardPairHMM::ForwardPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dict, Definitions::ModelType model,
		bool banding, unsigned int bandPercentage, unsigned int rateCategories, Maths* mth, Definitions::DpMatrixType mt) :
		EvolutionaryPairHMM(s1,s2, dict, rateCategories, mth, model, banding, bandPercentage,mt)
{
}

ForwardPairHMM::~ForwardPairHMM()
{
}

double ForwardPairHMM::runAlgorithm()
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

	for (i = 0; i<xSize; i++)
	{
		for (j = 0; j<ySize; j++)
		{
			if(this->withinBand(i,j,this->bandSpan) || !bandingEnabled)
			{
				if(i!=0)
				{
					k = i-1;
					emissionX = log(ptmatrix->getEquilibriumFreq(seq1[i-1].getMatrixIndex()));
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
	}

	sM = M->getValueAt(xSize-1, ySize-1);
	sX = X->getValueAt(xSize-1, ySize-1);
	sY = Y->getValueAt(xSize-1, ySize-1);
	sS = maths->logSum(sM,sX,sY);

	//cerr << "\t" << sX << "\t" << sY << "\t"<< sM << "\t" << sS << endl;

	//DEBUG ("Forward results:");
	//DEBUG (" sX, sY, sM, sS " << sX << "\t" << sY << "\t" << sM << "\t" << sS);

	return sS* -1.0;
}


} /* namespace EBC */
