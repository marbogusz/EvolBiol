/*
 * EvolutionaryPairHMM.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#include "hmm/EvolutionaryPairHMM.hpp"
#include "hmm/DpMatrixFull.hpp"
#include "models/NegativeBinomialGapModel.hpp"

namespace EBC
{

EvolutionaryPairHMM::EvolutionaryPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2,
			SubstitutionModelBase* smdl, IndelModel* imdl,
			Definitions::DpMatrixType mt, Band* bandObj) : substModel(smdl), indelModel(imdl), band(bandObj)
{
	M = X = Y = NULL;

	this->seq1 = s1;
	this->seq2 = s2;

	this->xSize = seq1.size() +1;
	this->ySize = seq2.size() +1;

	DUMP("Evolutionary Pair HMM for seqence 1 with size " << xSize << " and sequence 2 with size " << ySize);

	ptmatrix = new PMatrixDouble(substModel);

	this->tpb = new TransitionProbabilities(indelModel);

	piM = 0;
	piI = Definitions::minMatrixLikelihood;
	piD = Definitions::minMatrixLikelihood;

	initializeStates(mt);
}

void EvolutionaryPairHMM::setDivergenceTime(double time)
{
	ptmatrix->setTime(time);
	tpb->setTime(time);


}

void EvolutionaryPairHMM::getStateEquilibriums()
{
	md[0][0] = 1.0-2*g;
	md[1][1] = e+((1.0-e)*g);
	md[2][2] = e+((1.0-e)*g);
	md[0][1] = g;
	md[0][2] = g;

	md[1][0] = (1.0-e)*(1-2*g);
	md[2][0] = (1.0-e)*(1-2*g);

	md[2][1] = (1.0-e)*g;
	md[1][2] = (1.0-e)*g;

	piD = ((1.0-md[0][0])+(md[0][1]*(1.0-md[0][0]+md[1][0])/(md[1][1]-1.0-md[0][1])))/(((md[0][1]-md[2][1])*(1.0-md[0][0]+md[1][0])/(md[1][1]-1.0-md[0][1]))+md[2][0]-md[0][0]+1);
	piI = ((piD*(md[0][1]-md[2][1]))-md[0][1])/(md[1][1]-1.0-md[0][1]);
	piM = 1.0 -piI - piD;
}

void EvolutionaryPairHMM::setTransitionProbabilities()
{
	e = tpb->getGapExtension();
	g = tpb->getGapOpening();

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

void EvolutionaryPairHMM::summarize()
{
	e = tpb->getGapExtension();
	g = tpb->getGapOpening();

	cout << " Transition probabilities: " << endl;
	cout << "M->M : " << 1-2*g << endl;
	cout << "I->I : " << e+((1-e)*g) << endl;
	cout << "M->I : " << g << endl;
	cout << "I->M : " << (1-2*g)*(1-e) << endl;
	cout << "I->D : " << (1-e)*g << endl << endl;

	indelModel->summarize();
	substModel->summarize();

}

void EvolutionaryPairHMM::initializeStates(Definitions::DpMatrixType mt)
{

	if (M != NULL)
		delete M;
	if (X != NULL)
		delete X;
	if (Y != NULL)
		delete Y;

	switch (mt)
	{
	case Definitions::DpMatrixType::Full :
		M = new PairwiseHmmMatchState(xSize,ySize);
		X = new PairwiseHmmInsertState(xSize,ySize);
		Y = new PairwiseHmmDeleteState(xSize,ySize);
		break;
	case Definitions::DpMatrixType::Limited :
		M = new PairwiseHmmMatchState(new DpMatrixLoMem(xSize,ySize));
		X = new PairwiseHmmInsertState(new DpMatrixLoMem(xSize,ySize));
		Y = new PairwiseHmmDeleteState(new DpMatrixLoMem(xSize,ySize));
		break;
	default :
		M = new PairwiseHmmMatchState(xSize,ySize);
		X = new PairwiseHmmInsertState(xSize,ySize);
		Y = new PairwiseHmmDeleteState(xSize,ySize);
	}
}

void EvolutionaryPairHMM::calculateModels()
{
	ptmatrix->calculate();
	tpb->calculate();
}

EvolutionaryPairHMM::~EvolutionaryPairHMM()
{

	DUMP("Evolutionary pair HMM destructor");
	delete Y;
	delete X;
	delete M;
    delete ptmatrix;
    delete tpb;
}

double EvolutionaryPairHMM::getAlignmentLikelihood(vector<SequenceElement> s1,
		vector<SequenceElement> s2)
{
	double lnl = 0;

	//cerr << endl;

	int k =0;
	int l =0;
	PairwiseHmmStateBase* previous;
	if(s2[0].isIsGap()){
		previous = X;
		k++;
		lnl += ptmatrix->getLogEquilibriumFreq(s1[0].getMatrixIndex());
		//DUMP("I " << 0 << "\tlnl\t" << lnl << "\tmatrix\t" << previous->getValueAt(k,l));
	}
	else if(s1[0].isIsGap()){
		previous = Y;
		l++;
		lnl += ptmatrix->getLogEquilibriumFreq(s2[0].getMatrixIndex());
		//DUMP("D " << 0 << "\tlnl\t" << lnl << "\tmatrix\t" << previous->getValueAt(k,l));
	}
	else{
		previous = M;
		k++;
		l++;
		lnl += ptmatrix->getLogPairTransition(s1[0].getMatrixIndex(), s2[0].getMatrixIndex());
		//DUMP("M " << 0 << "\tlnl\t" << lnl << "\tmatrix\t" << previous->getValueAt(k,l));
	}
	for(int i=1; i< s1.size(); i++){
		if(s2[i].isIsGap()){
			//Insert
			k++;
			if (previous == X)
				lnl += X->getTransitionProbabilityFromInsert();
			else if (previous == Y)
				lnl += X->getTransitionProbabilityFromDelete();
			else
				lnl+=X->getTransitionProbabilityFromMatch();
			lnl += ptmatrix->getLogEquilibriumFreq(s1[i].getMatrixIndex());
			//DUMP("I " 0 << "\tlnl\t" << lnl << "\tmatrix\t" << previous->getValueAt(i,j));
			previous = X;
		}
		else if(s1[i].isIsGap()){
			//Delete
			l++;
			if (previous == X)
				lnl += Y->getTransitionProbabilityFromInsert();
			else if (previous == Y)
				lnl += Y->getTransitionProbabilityFromDelete();
			else
				lnl+=Y->getTransitionProbabilityFromMatch();
			lnl += ptmatrix->getLogEquilibriumFreq(s2[i].getMatrixIndex());
			//DUMP("D " 0 << "\tlnl\t" << lnl << "\tmatrix\t" << previous->getValueAt(i,j));
			previous = Y;
		}
		else{
			//Match
			k++;
			l++;
			if (previous == X)
				lnl += M->getTransitionProbabilityFromInsert();
			else if (previous == Y)
				lnl += M->getTransitionProbabilityFromDelete();
			else
				lnl+=M->getTransitionProbabilityFromMatch();
			lnl += ptmatrix->getLogPairTransition(s1[i].getMatrixIndex(), s2[i].getMatrixIndex());
			//DUMP("M " 0 << "\tlnl\t" << lnl << "\tmatrix\t" << previous->getValueAt(i,j));
			previous = M;
		}
		//DUMP(i << "\tlnl " << lnl);
	}
	//cerr << endl;
	return lnl;

}

} /* namespace EBC */


