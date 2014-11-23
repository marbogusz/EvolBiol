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


ViterbiPairHMM::ViterbiPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, SubstitutionModelBase* smdl,
		IndelModel* imdl,Definitions::DpMatrixType mt, Band* bandObj) :
		EvolutionaryPairHMM(s1,s2, smdl, imdl,mt, bandObj)
{
	this->alignment.reserve(xSize);
}

ViterbiPairHMM::~ViterbiPairHMM()
{
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

double ViterbiPairHMM::runAlgorithm()
{

	calculateModels();
	setTransitionProbabilities();

	unsigned int i,j,k,l;

	double xx,xy,xm,yx,yy,ym,mx,my,mm;

	double emissionM;
	double emissionX;
	double emissionY;

	M->initializeData();
	X->initializeData();
	Y->initializeData();

		//while (i != xSize && j != ySize)

	for (i = 0; i<xSize; i++)
	{
		for (j = 0; j<ySize; j++)
		{
			if(i!=0)
			{

				k = i-1;
				emissionX = ptmatrix->getLogEquilibriumFreq(seq1[i-1].getMatrixIndex());
				xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
				xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
				xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();

				X->setValueAt(i,j,getMax(xm,xx,xy,i,j,X) + emissionX);
			}
			if(j!=0)
			{
				k = j-1;
				emissionY = ptmatrix->getLogEquilibriumFreq(seq2[j-1].getMatrixIndex());
				ym = M->getValueAt(i,k) + Y->getTransitionProbabilityFromMatch();
				yx = X->getValueAt(i,k) + Y->getTransitionProbabilityFromInsert();
				yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
				Y->setValueAt(i,j,getMax(ym,yx,yy,i,j,Y) + emissionY);
			}

			if(i!=0 && j!=0)
			{
				k = i-1;
				l = j-1;
				emissionM = ptmatrix->getLogPairTransition(seq1[i-1].getMatrixIndex(), seq2[j-1].getMatrixIndex());
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
		M->tracebackRaw(this->seq1,this->seq2, this->dict, this->alignment);
	}
	else if(mx >= my)
	{
		X->tracebackRaw(this->seq1,this->seq2, this->dict, this->alignment);
	}
	else
	{
		Y->tracebackRaw(this->seq1,this->seq2, this->dict, this->alignment);
	}
*/
	DUMP("Final Viterbi M  " << mm);
	DUMP("Final Viterbi X  " << mx );
	DUMP("Final Viterbi Y  " << my );

	this->setTotalLikelihood(std::max(mm,std::max(mx,my)));

return getTotalLikelihood()*-1.0;
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
