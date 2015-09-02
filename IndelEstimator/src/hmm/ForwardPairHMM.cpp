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


ForwardPairHMM::ForwardPairHMM(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, SubstitutionModelBase* smdl,
		IndelModel* imdl, Definitions::DpMatrixType mt, Band* bandObj, bool useEquilibriumFreqs) :
		EvolutionaryPairHMM(s1,s2, smdl, imdl, mt, bandObj, true)
{
}

ForwardPairHMM::~ForwardPairHMM()
{
}


pair<string, string> ForwardPairHMM::getBestAlignment(string&seq_a, string& seq_b)
{
	DUMP("Forward HMM getBestAlignment");
	pair<string, string> alignment;

	//reserve memory for out strings (20% of gaps should be ok)
	alignment.first.reserve(max(xSize,ySize)*1.2);
	alignment.second.reserve(max(xSize,ySize)*1.2);


	//std::random_device rd;
	//std::mt19937 gen(rd());
	//std::uniform_real_distribution<> dis(0, 1.0);

	unsigned int i = xSize-1;
	unsigned int j = ySize-1;

	double mtProb,inProb,dlProb,currProb;
	double emission = 0.0;

	//choose initial state
	PairwiseHmmStateBase* currentState;

	mtProb = M->getValueAt(xSize-1, ySize-1) - this->getTotalLikelihood();
	inProb = X->getValueAt(xSize-1, ySize-1) - this->getTotalLikelihood();
	dlProb = Y->getValueAt(xSize-1, ySize-1) - this->getTotalLikelihood();

	while(i > 0 && j > 0)
	{

		if(mtProb >= inProb && mtProb >= dlProb )
			currentState = M;
		else if (inProb >= dlProb)
			currentState = X;
		else currentState = Y;

		currProb = currentState->getValueAt(i,j);
		if (currentState->stateId == Definitions::StateId::Match)
		{
			emission = ptmatrix->getLogPairTransition((*seq1)[i-1]->getMatrixIndex(), (*seq2)[j-1]->getMatrixIndex());
			alignment.first += seq_a[i-1];
			alignment.second += seq_b[j-1];
			i--;
			j--;
		}
		else if (currentState->stateId == Definitions::StateId::Delete)
		{
			ptmatrix->getLogEquilibriumFreq((*seq2)[j-1]->getMatrixIndex());
			alignment.second += seq_b[j-1];
			alignment.first += '-';
			j--;
		}
		else //Insert
		{
			emission = ptmatrix->getLogEquilibriumFreq((*seq1)[i-1]->getMatrixIndex());
			alignment.first += seq_a[i-1];
			alignment.second += '-';
			i--;
		}
		mtProb = emission + currentState->getTransitionProbabilityFromMatch() + M->getValueAt(i,j) - currProb;
		inProb = emission + currentState->getTransitionProbabilityFromInsert() + X->getValueAt(i,j) - currProb;
		dlProb = emission + currentState->getTransitionProbabilityFromDelete() + Y->getValueAt(i,j) - currProb;
	}

	if (j==0)
	{
		while(i > 0){
			alignment.first += seq_a[i-1];
			alignment.second += '-';
			i--;
		}
	}
	else if (i==0)
	{
		while(j > 0){
			alignment.second += seq_b[j-1];
			alignment.first += '-';
			j--;
		}
	}
	//deal with the last row or column

	reverse(alignment.first.begin(), alignment.first.end());
	reverse(alignment.second.begin(), alignment.second.end());
	return alignment;
}

double ForwardPairHMM::runAlgorithm()
{

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

	//DUMP("Forward equilibriums : PiM\t" << piM << "\tPiI\t" << piI << "\tPiD\t" << piD);

	M->initializeData(this->piM);
	X->initializeData(this->piI);
	Y->initializeData(this->piD);

	if(this->band == NULL)
	{
		//handle 0-index rows and columns separately!
		//1st col
		X->setValueAt(1,0, ptmatrix->getLogEquilibriumFreq((*seq1)[0]->getMatrixIndex()) + initTransX);

		for(i=2,j=0; i< xSize; i++)
		{
			k = i-1;
			emissionX = ptmatrix->getLogEquilibriumFreq((*seq1)[i-1]->getMatrixIndex());
			xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
			X->setValueAt(i,j, emissionX + xx);
		}
		//1st row

		Y->setValueAt(0,1, ptmatrix->getLogEquilibriumFreq((*seq2)[0]->getMatrixIndex()) + initTransY);
		for(j=2,i=0; j< ySize; j++)
		{
			k = j-1;
			emissionY = ptmatrix->getLogEquilibriumFreq((*seq2)[j-1]->getMatrixIndex());
			yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
			Y->setValueAt(i,j, emissionY + yy);
		}


		for (i = 1; i<xSize; i++)
		{
			for (j = 1; j<ySize; j++)
			{

				k = i-1;
				emissionX = ptmatrix->getLogEquilibriumFreq((*seq1)[i-1]->getMatrixIndex());
				xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
				xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
				xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();
				X->setValueAt(i,j, emissionX + maths->logSum(xm,xx,xy));

				k = j-1;
				emissionY = ptmatrix->getLogEquilibriumFreq((*seq2)[j-1]->getMatrixIndex());
				ym = M->getValueAt(i,k) + Y->getTransitionProbabilityFromMatch();
				yx = X->getValueAt(i,k) + Y->getTransitionProbabilityFromInsert();
				yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
				Y->setValueAt(i,j, emissionY + maths->logSum(ym,yx,yy));

				k = i-1;
				l = j-1;
				emissionM = ptmatrix->getLogPairTransition((*seq1)[i-1]->getMatrixIndex(), (*seq2)[j-1]->getMatrixIndex());
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
				emissionX = ptmatrix->getLogEquilibriumFreq((*seq1)[i-1]->getMatrixIndex());
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

			if (loD > -1)
			{
				hiD = bracketD.second;
				for(i = loD; i <= hiD; i++)
				{
					k = j-1;
					emissionY = ptmatrix->getLogEquilibriumFreq((*seq2)[j-1]->getMatrixIndex());
					ym = M->getValueAt(i,k) + Y->getTransitionProbabilityFromMatch();
					yx = X->getValueAt(i,k) + Y->getTransitionProbabilityFromInsert();
					yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
					Y->setValueAt(i,j, emissionY + maths->logSum(ym,yx,yy));
				}
			}
			if (loM > 0)
			{
				hiM = bracketM.second;
				for(i = loM; i <= hiM; i++)
				{
					k = i-1;
					l = j-1;
					emissionM = ptmatrix->getLogPairTransition((*seq1)[i-1]->getMatrixIndex(), (*seq2)[j-1]->getMatrixIndex());
					mm = M->getValueAt(k,l) + M->getTransitionProbabilityFromMatch();
					mx = X->getValueAt(k,l) + M->getTransitionProbabilityFromInsert();
					my = Y->getValueAt(k,l) + M->getTransitionProbabilityFromDelete();
					M->setValueAt(i,j, emissionM + maths->logSum(mm,mx,my));
				}
			}

			if (loI > 0)
			{
				hiI = bracketI.second;
				for(i = loI; i <= hiI; i++)
				{
					k = i-1;
					emissionX = ptmatrix->getLogEquilibriumFreq((*seq1)[i-1]->getMatrixIndex());
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

	sS = maths->logSum(sM,sX,sY) + log(xi);

	this->setTotalLikelihood(sS);

	//cerr << "\t" << sX << "\t" << sY << "\t"<< sM << "\t" << sS << endl;

	DUMP ("Forward lnls I, D, M, Total " << sX << "\t" << sY << "\t" << sM << "\t" << sS);
	DEBUG(sS << " total forward likelihood");

	return sS* -1.0;
}



} /* namespace EBC */
