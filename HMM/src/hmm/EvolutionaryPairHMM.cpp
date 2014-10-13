/*
 * EvolutionaryPairHMM.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#include "hmm/EvolutionaryPairHMM.hpp"
#include "models/NegativeBinomialGapModel.hpp"

namespace EBC
{

EvolutionaryPairHMM::EvolutionaryPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, bool banding,
			SubstitutionModelBase* smdl, IndelModel* imdl, unsigned int bandPercentage, Definitions::DpMatrixType mt) :
				substModel(smdl), indelModel(imdl)
{
	M = X = Y = NULL;

	this->seq1 = s1;
	this->seq2 = s2;

	this->xSize = seq1.size() +1;
	this->ySize = seq2.size() +1;

	ptmatrix = new PMatrixDouble(substModel);

	bandFactor = bandPercentage;
	bandingEnabled = banding;

	this->tpb = new TransitionProbabilities(indelModel);

	getBandWidth();
	initializeStates(mt);

}

void EvolutionaryPairHMM::setDivergenceTime(double time)
{
	ptmatrix->setTime(time);
	tpb->setTime(time);


}

void EvolutionaryPairHMM::setTransitionProbabilities()
{
	double e,g;

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
	double e,g;
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

	delete Y;
	delete X;
	delete M;
    delete ptmatrix;
}

} /* namespace EBC */
