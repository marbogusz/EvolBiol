/*
 * ViterbiPairHMM.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "Definitions.hpp"
#include "ViterbiPairHMM.hpp"
#include <algorithm>

namespace EBC
{


ViterbiPairHMM::ViterbiPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dict, Definitions::ModelType model,
		bool banding, unsigned int bandPercentage, unsigned int rateCategories, Maths* mt) :
		EvolutionaryPairHMM(s1,s2, dict, rateCategories, mt, model, banding, bandPercentage)
{
}

ViterbiPairHMM::~ViterbiPairHMM()
{
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

double ViterbiPairHMM::runViterbiAlgorithm()
{

	calculateModels();
	setTransitionProbabilities();

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
	mx = X->getValueAt(xSize-1,ySize-1);
	my = Y->getValueAt(xSize-1,ySize-1);
	mm = M->getValueAt(xSize-1,ySize-1);

	DEBUG("Final Viterbi M  " << mm);
	DEBUG("Final Viterbi X  " << mx );
	DEBUG("Final Viterbi Y  " << my );

return std::max(mm,std::max(mx,my));

}


} /* namespace EBC */
