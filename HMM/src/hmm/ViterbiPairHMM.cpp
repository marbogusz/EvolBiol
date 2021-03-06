//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================

/*
 * ViterbiPairHMM.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "core/Definitions.hpp"
#include "hmm/ViterbiPairHMM.hpp"
#include <algorithm>

namespace EBC
{


ViterbiPairHMM::ViterbiPairHMM(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, SubstitutionModelBase* smdl,
		IndelModel* imdl,Definitions::DpMatrixType mt, Band* bandObj, bool useEquilibriumFreqs ) :
		EvolutionaryPairHMM(s1,s2, smdl, imdl,mt, bandObj, useEquilibriumFreqs)
{
	this->alignment.reserve(xSize);
}

ViterbiPairHMM::~ViterbiPairHMM()
{
}

void ViterbiPairHMM::getAlignment(HMMPathSample& sample)
{

	unsigned char gapElem = this->substModel->getMatrixSize(); // last matrix element is the gap ID!

	unsigned int i = xSize-1;
	unsigned int j = ySize-1;

	double mtProb, inProb, dlProb, currProb, rnbr, tmp;
	double emission = 0.0;

	//choose initial state
	PairwiseHmmStateBase* currentState;
	PairwiseHmmStateBase* previousState;

	mtProb = M->getValueAt(xSize-1, ySize-1);
	inProb = X->getValueAt(xSize-1, ySize-1);
	dlProb = Y->getValueAt(xSize-1, ySize-1);



	//Establish which state is first!
	//TODO - remove the mess after testing!

	previousState = NULL;

	if(mtProb >= inProb && mtProb >= dlProb )
		currentState = M;
	else if (inProb >= dlProb)
		currentState = X;
	else currentState = Y;

	currProb = currentState->getValueAt(i,j);

	if (currentState->stateId == Definitions::StateId::Match)
	{
		emission = ptmatrix->getLogPairTransition((*seq1)[i-1]->getMatrixIndex(), (*seq2)[j-1]->getMatrixIndex());
		sample.addSitePattern((*seq1)[i-1]->getMatrixIndex(),(*seq2)[j-1]->getMatrixIndex());
		i--;
		j--;
	}
	else if (currentState->stateId == Definitions::StateId::Delete)
	{
		ptmatrix->getLogEquilibriumFreq((*seq2)[j-1]->getMatrixIndex());
		sample.addSitePattern(gapElem,(*seq2)[j-1]->getMatrixIndex());
		j--;
	}
	else //Insert
	{
		emission = ptmatrix->getLogEquilibriumFreq((*seq1)[i-1]->getMatrixIndex());
		sample.addSitePattern((*seq1)[i-1]->getMatrixIndex(),gapElem);
		i--;
	}
	mtProb = currentState->getTransitionProbabilityFromMatch() + M->getValueAt(i,j);
	inProb = currentState->getTransitionProbabilityFromInsert() + X->getValueAt(i,j);
	dlProb = currentState->getTransitionProbabilityFromDelete() + Y->getValueAt(i,j);

	previousState = currentState;

	//Calculate the rest;

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
			sample.addSitePattern((*seq1)[i-1]->getMatrixIndex(),(*seq2)[j-1]->getMatrixIndex());
			i--;
			j--;
		}
		else if (currentState->stateId == Definitions::StateId::Delete)
		{
			ptmatrix->getLogEquilibriumFreq((*seq2)[j-1]->getMatrixIndex());
			sample.addSitePattern(gapElem,(*seq2)[j-1]->getMatrixIndex());
			j--;
		}
		else //Insert
		{
			emission = ptmatrix->getLogEquilibriumFreq((*seq1)[i-1]->getMatrixIndex());
			sample.addSitePattern((*seq1)[i-1]->getMatrixIndex(),gapElem);
			i--;
		}
		mtProb = currentState->getTransitionProbabilityFromMatch() + M->getValueAt(i,j);
		inProb = currentState->getTransitionProbabilityFromInsert() + X->getValueAt(i,j);
		dlProb = currentState->getTransitionProbabilityFromDelete() + Y->getValueAt(i,j);

		sample.addTransition(currentState->stateId, previousState->stateId );

		previousState = currentState;

	}

	if (j==0)
	{
		while(i > 0){
			currentState = X;
			sample.addSitePattern((*seq1)[i-1]->getMatrixIndex(),gapElem);
			sample.addTransition(currentState->stateId, previousState->stateId );
			previousState = currentState;
			i--;
		}
	}
	else if (i==0)
	{
		while(j > 0){
			currentState = Y;
			sample.addSitePattern(gapElem,(*seq2)[j-1]->getMatrixIndex());
			sample.addTransition(currentState->stateId, previousState->stateId );
			previousState = currentState;
			j--;
		}
	}
	sample.setLastState(currentState);
}

pair<string, string> ViterbiPairHMM::getBestAlignment(string&seq_a, string& seq_b)
{
	DUMP("Viterbi HMM get alignment");
	pair<string, string> alignment;

	//reserve memory for out strings (20% of gaps should be ok)
	alignment.first.reserve(max(xSize,ySize)*1.2);
	alignment.second.reserve(max(xSize,ySize)*1.2);


	unsigned int i = xSize-1;
	unsigned int j = ySize-1;

	double mtProb,inProb,dlProb,currProb, rnbr,tmp;
	double emission;

	//choose initial state
	PairwiseHmmStateBase* currentState;

	mtProb = M->getValueAt(xSize-1, ySize-1);
	inProb = X->getValueAt(xSize-1, ySize-1);
	dlProb = Y->getValueAt(xSize-1, ySize-1);

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
		mtProb = currentState->getTransitionProbabilityFromMatch() + M->getValueAt(i,j);
		inProb = currentState->getTransitionProbabilityFromInsert() + X->getValueAt(i,j);
		dlProb = currentState->getTransitionProbabilityFromDelete() + Y->getValueAt(i,j);
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

pair<string, string> ViterbiPairHMM::getStringAlignment()
{
	DUMP("Viterbi HMM get string alignment");
	pair<string, string> alignment;

	//reserve memory for out strings (20% of gaps should be ok)
	alignment.first.reserve(max(xSize,ySize)*1.2);
	alignment.second.reserve(max(xSize,ySize)*1.2);


	unsigned int i = xSize-1;
	unsigned int j = ySize-1;

	double mtProb,inProb,dlProb,currProb, rnbr,tmp;
	double emission;

	//choose initial state
	PairwiseHmmStateBase* currentState;

	mtProb = M->getValueAt(xSize-1, ySize-1);
	inProb = X->getValueAt(xSize-1, ySize-1);
	dlProb = Y->getValueAt(xSize-1, ySize-1);

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
			alignment.first += (*seq1)[i-1]->getSymbol();
			alignment.second += (*seq2)[j-1]->getSymbol();
			i--;
			j--;
		}
		else if (currentState->stateId == Definitions::StateId::Delete)
		{
			ptmatrix->getLogEquilibriumFreq((*seq2)[j-1]->getMatrixIndex());
			alignment.second += (*seq2)[j-1]->getSymbol();
			alignment.first += '-';
			j--;
		}
		else //Insert
		{
			emission = ptmatrix->getLogEquilibriumFreq((*seq1)[i-1]->getMatrixIndex());
			alignment.first += (*seq1)[i-1]->getSymbol();
			alignment.second += '-';
			i--;
		}
		mtProb = currentState->getTransitionProbabilityFromMatch() + M->getValueAt(i,j);
		inProb = currentState->getTransitionProbabilityFromInsert() + X->getValueAt(i,j);
		dlProb = currentState->getTransitionProbabilityFromDelete() + Y->getValueAt(i,j);
	}

	if (j==0)
	{
		while(i > 0){
			alignment.first += (*seq1)[i-1]->getSymbol();
			alignment.second += '-';
			i--;
		}
	}
	else if (i==0)
	{
		while(j > 0){
			alignment.second += (*seq2)[j-1]->getSymbol();
			alignment.first += '-';
			j--;
		}
	}
	//deal with the last row or column

	reverse(alignment.first.begin(), alignment.first.end());
	reverse(alignment.second.begin(), alignment.second.end());
	return alignment;
}

double ViterbiPairHMM::getViterbiSubstitutionLikelihood()
{
	double lnl = 0;
	calculateModels();
	//we have the site patterns now
	//get alignment!

	for (auto it = alignment.begin(); it != alignment.end(); it ++)
	{
		lnl += this->ptmatrix->getPairSitePattern(it->first, it->second);
	}

	return lnl * -1.0;
}

double ViterbiPairHMM::getMax(double m, double x, double y, unsigned int i, unsigned int j, PairwiseHmmStateBase* state)
{
	if(m >x && m >y)
	{
		//state->setDiagonalAt(i,j);

		//cout << i << " " << j << " coming from M" << endl;
		state->setDirection(i,j);
		state->setSourceMatrixPtr(i,j,M);
		return m;
	}
	else if(x > y)
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

double ViterbiPairHMM::runAlgorithm()
{
	unsigned int i,j,k,l;

	double xx,xy,xm,yx,yy,ym,mx,my,mm, sS;

	double emissionM;
	double emissionX;
	double emissionY;

	DUMP("Viterbi equilibriums : PiM\t" << piM << "\tPiI\t" << piI << "\tPiD\t" << piD);

	M->initializeData(this->piM);
	X->initializeData(this->piI);
	Y->initializeData(this->piD);

		//while (i != xSize && j != ySize)

	for (i = 0; i<xSize; i++)
	{
		for (j = 0; j<ySize; j++)
		{
			if(i!=0)
			{

				k = i-1;
				emissionX = ptmatrix->getLogEquilibriumFreq((*seq1)[i-1]->getMatrixIndex());
				xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
				xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
				xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();

				X->setValueAt(i,j,getMax(xm,xx,xy,i,j,X) + emissionX);
			}
			if(j!=0)
			{
				k = j-1;
				emissionY = ptmatrix->getLogEquilibriumFreq((*seq2)[j-1]->getMatrixIndex());
				ym = M->getValueAt(i,k) + Y->getTransitionProbabilityFromMatch();
				yx = X->getValueAt(i,k) + Y->getTransitionProbabilityFromInsert();
				yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
				Y->setValueAt(i,j,getMax(ym,yx,yy,i,j,Y) + emissionY);
			}

			if(i!=0 && j!=0)
			{
				k = i-1;
				l = j-1;
				emissionM = ptmatrix->getLogPairTransition((*seq1)[i-1]->getMatrixIndex(), (*seq2)[j-1]->getMatrixIndex());
				mm = M->getValueAt(k,l) + M->getTransitionProbabilityFromMatch();
				mx = X->getValueAt(k,l) + M->getTransitionProbabilityFromInsert();
				my = Y->getValueAt(k,l) + M->getTransitionProbabilityFromDelete();
				M->setValueAt(i,j,getMax(mm,mx,my,i,j,M) + emissionM);
			}
		}
	}
	mx = X->getValueAt(xSize-1,ySize-1);
	my = Y->getValueAt(xSize-1,ySize-1);
	mm = M->getValueAt(xSize-1,ySize-1);

/*
	if(mm >=mx && mm >=my)
	{
		M->tracebackRaw(this->(*seq1),this->(*seq2), this->dict, this->alignment);
	}
	else if(mx >= my)
	{
		X->tracebackRaw(this->(*seq1),this->(*seq2), this->dict, this->alignment);
	}
	else
	{
		Y->tracebackRaw(this->(*seq1),this->(*seq2), this->dict, this->alignment);
	}
*/
	DUMP("Final Viterbi M  " << mm);
	DUMP("Final Viterbi X  " << mx );
	DUMP("Final Viterbi Y  " << my );

return (std::max(mm,std::max(mx,my)))*-1.0;
}

pair<string, string> ViterbiPairHMM::getAlignment(string&a, string& b)
{
	double mv, xv, yv;

	pair<string, string> initialAlignment;
	initialAlignment.first.reserve(xSize > ySize ? xSize : ySize);
	initialAlignment.second.reserve(xSize > ySize ? xSize : ySize);




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


	return initialAlignment;
}


} /* namespace EBC */
