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


ForwardPairHMM::ForwardPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, SubstitutionModelBase* smdl,
		IndelModel* imdl, Definitions::DpMatrixType mt, Band* bandObj) :
		EvolutionaryPairHMM(s1,s2, smdl, imdl, mt, bandObj)
{
}

ForwardPairHMM::~ForwardPairHMM()
{
}

double ForwardPairHMM::runAlgorithm()
{

	calculateModels();
	setTransitionProbabilities();

	int i;
	int j;
	int k;
	int l;

	double sX,sY,sM, sS;

	double xx,xy,xm,yx,yy,ym,mx,my,mm;

	double emissionM;
	double emissionX;
	double emissionY;

	//FIXME - multiple runs using the same hmm object do not require dp matrix zeroing as long as the band stays the same!
	M->initializeData();
	X->initializeData();
	Y->initializeData();

	if(this->band == NULL)
	{
		//handle 0-index rows and columns separately!
		//1st col
		for(i=1,j=0; i< xSize; i++)
		{
			k = i-1;
			emissionX = ptmatrix->getLogEquilibriumFreq(seq1[i-1].getMatrixIndex());
			xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
			xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
			xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();
			X->setValueAt(i,j, emissionX + maths->logSum(xm,xx,xy));
		}
		//1st row
		for(j=1,i=0; j< ySize; j++)
		{
			k = j-1;
			emissionY = ptmatrix->getLogEquilibriumFreq(seq2[j-1].getMatrixIndex());
			ym = M->getValueAt(i,k) + Y->getTransitionProbabilityFromMatch();
			yx = X->getValueAt(i,k) + Y->getTransitionProbabilityFromInsert();
			yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
			Y->setValueAt(i,j, emissionY + maths->logSum(ym,yx,yy));
		}


		for (i = 1; i<xSize; i++)
		{
			for (j = 1; j<ySize; j++)
			{

				k = i-1;
				emissionX = ptmatrix->getLogEquilibriumFreq(seq1[i-1].getMatrixIndex());
				xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
				xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
				xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();
				X->setValueAt(i,j, emissionX + maths->logSum(xm,xx,xy));

				k = j-1;
				emissionY = ptmatrix->getLogEquilibriumFreq(seq2[j-1].getMatrixIndex());
				ym = M->getValueAt(i,k) + Y->getTransitionProbabilityFromMatch();
				yx = X->getValueAt(i,k) + Y->getTransitionProbabilityFromInsert();
				yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
				Y->setValueAt(i,j, emissionY + maths->logSum(ym,yx,yy));

				k = i-1;
				l = j-1;
				emissionM = ptmatrix->getLogPairTransition(seq1[i-1].getMatrixIndex(), seq2[j-1].getMatrixIndex());
				mm = M->getValueAt(k,l) + M->getTransitionProbabilityFromMatch();
				mx = X->getValueAt(k,l) + M->getTransitionProbabilityFromInsert();
				my = Y->getValueAt(k,l) + M->getTransitionProbabilityFromDelete();
				M->setValueAt(i,j, emissionM + maths->logSum(mm,mx,my));
			}
		}
	}
	else
	{
		//banding column by column!
		//calculate 1st column for X
		int loI, hiI, loD, hiD, loM, hiM;
		auto bracket = band->getInsertRangeAt(0);
		loI = bracket.first;
		if (loI > 0)
		{
			hiI = bracket.second;
			for(i=loI,j=0; i<= hiI; i++)
			{
				k = i-1;
				emissionX = ptmatrix->getLogEquilibriumFreq(seq1[i-1].getMatrixIndex());
				xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
				xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
				xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();
				X->setValueAt(i,j, emissionX + maths->logSum(xm,xx,xy));
			}
		}
		for(j=1; j<ySize; j++)
		{
			//FIXME - range should be a reference perhaps
			auto bracketI = band->getInsertRangeAt(j);
			auto bracketD = band->getDeleteRangeAt(j);
			auto bracketM = band->getMatchRangeAt(j);
			loI = bracketI.first;
			loM = bracketM.first;
			loD = bracketD.first;

			if (loD != -1)
			{
				hiD = bracketD.second;
				for(i = loD; i <= hiD; i++)
				{
					k = j-1;
					emissionY = ptmatrix->getLogEquilibriumFreq(seq2[j-1].getMatrixIndex());
					ym = M->getValueAt(i,k) + Y->getTransitionProbabilityFromMatch();
					yx = X->getValueAt(i,k) + Y->getTransitionProbabilityFromInsert();
					yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
					Y->setValueAt(i,j, emissionY + maths->logSum(ym,yx,yy));
				}
			}
			if (loM != -1)
			{
				hiM = bracketM.second;
				for(i = loM == 0 ? 1 : loM; i <= hiM; i++)
				{
					k = i-1;
					l = j-1;
					emissionM = ptmatrix->getLogPairTransition(seq1[i-1].getMatrixIndex(), seq2[j-1].getMatrixIndex());
					mm = M->getValueAt(k,l) + M->getTransitionProbabilityFromMatch();
					mx = X->getValueAt(k,l) + M->getTransitionProbabilityFromInsert();
					my = Y->getValueAt(k,l) + M->getTransitionProbabilityFromDelete();
					M->setValueAt(i,j, emissionM + maths->logSum(mm,mx,my));
				}
			}

			if (loI != -1)
			{
				hiI = bracketI.second;
				for(i = loI == 0 ? 1 : loI; i <= hiI; i++)
				{
					k = i-1;
					emissionX = ptmatrix->getLogEquilibriumFreq(seq1[i-1].getMatrixIndex());
					xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
					xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
					xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();
					X->setValueAt(i,j, emissionX + maths->logSum(xm,xx,xy));
				}
			}
		}
	}

	sM = M->getValueAt(xSize-1, ySize-1);
	sX = X->getValueAt(xSize-1, ySize-1);
	sY = Y->getValueAt(xSize-1, ySize-1);
	sS = maths->logSum(sM,sX,sY);

	//cerr << "\t" << sX << "\t" << sY << "\t"<< sM << "\t" << sS << endl;

	DUMP ("Forward results:");
	DUMP (" sX, sY, sM, sS " << sX << "\t" << sY << "\t" << sM << "\t" << sS);

	return sS* -1.0;
}


} /* namespace EBC */
