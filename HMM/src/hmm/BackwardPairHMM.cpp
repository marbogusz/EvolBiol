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


BackwardPairHMM::BackwardPairHMM(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, SubstitutionModelBase* smdl,
		IndelModel* imdl,  Definitions::DpMatrixType mt ,Band* bandObj) :
		EvolutionaryPairHMM(s1,s2, smdl, imdl, mt, bandObj, true)
{
	this->MPstate = NULL;
}

BackwardPairHMM::~BackwardPairHMM()
{
	if (MPstate != NULL)
		delete MPstate;
}

void BackwardPairHMM::calculatePosteriors(ForwardPairHMM* fwd)
{
	DEBUG("Calculating posterior probabilities");

	unsigned int i,j;
	double xval, yval, mval;
	//totalForward
	double fwdT;

	fwdT = fwd->getTotalLikelihood();

/*
	DUMP("MATCH");
	//DUMP("FORWARD MATRICES");
	//dynamic_cast<DpMatrixFull*>(fwd->M->getDpMatrix())->outputValues(0);
	DUMP("BACKWARD MATRICES");
	dynamic_cast<DpMatrixFull*>(M->getDpMatrix())->outputValues(0);
	DUMP("\nINSERT");
	//DUMP("FORWARD MATRICES");
	//dynamic_cast<DpMatrixFull*>(fwd->X->getDpMatrix())->outputValues(0);
	DUMP("BACKWARD MATRICES");
	dynamic_cast<DpMatrixFull*>(X->getDpMatrix())->outputValues(0);

	DUMP("\nDELETE");
	//DUMP("FORWARD MATRICES");
	//dynamic_cast<DpMatrixFull*>(fwd->Y->getDpMatrix())->outputValues(0);
	DUMP("BACKWARD MATRICES");
	dynamic_cast<DpMatrixFull*>(Y->getDpMatrix())->outputValues(0);
	//DUMP("#####Match posteriors########");

	//DUMP("#####Insert posteriors########");

	//dynamic_cast<DpMatrixFull*>(X->getDpMatrix())->outputValues(0);
	//DUMP("#####Delete posteriors########");

	//dynamic_cast<DpMatrixFull*>(Y->getDpMatrix())->outputValues(0);

	DUMP("BACKWARD MATRICES");
*/
	for (i = 0; i<xSize-1; i++)
	{
		for (j = 0; j<ySize-1; j++)
		{
			xval = X->getValueAt(i,j) + fwd->X->getValueAt(i,j) - fwdT;
			yval = Y->getValueAt(i,j) + fwd->Y->getValueAt(i,j) - fwdT;
			mval = M->getValueAt(i,j) + fwd->M->getValueAt(i,j) - fwdT;

			X->setValueAt(i,j,xval);
			Y->setValueAt(i,j,yval);
			M->setValueAt(i,j,mval);
		}
	}

	DUMP("#####Match posteriors########");
	dynamic_cast<DpMatrixFull*>(M->getDpMatrix())->outputValues(0);
	DUMP("#####Insert posteriors########");
	dynamic_cast<DpMatrixFull*>(X->getDpMatrix())->outputValues(0);
	DUMP("#####Delete posteriors########");
	dynamic_cast<DpMatrixFull*>(Y->getDpMatrix())->outputValues(0);

}


double BackwardPairHMM::getAlignmentLikelihood(vector<SequenceElement*>* s1,
		vector<SequenceElement*>* s2, bool post, vector<vector<double> >& posteriors)
{
	double lnl = 0;
	double avp = 0;
	int k =0;
	int l =0;
	PairwiseHmmStateBase* previous;
	if((*s2)[0]->isIsGap()){
		previous = X;
		k++;
		lnl += ptmatrix->getLogEquilibriumFreq((*s1)[0]->getMatrixIndex()) + initTransX;
	}
	else if((*s1)[0]->isIsGap()){
		previous = Y;
		l++;
		lnl += ptmatrix->getLogEquilibriumFreq((*s2)[0]->getMatrixIndex()) + initTransY;
	}
	else{
		previous = M;
		k++;
		l++;
		posteriors[k][l] += 1.0;
		lnl += ptmatrix->getLogPairTransition((*s1)[0]->getMatrixIndex(), (*s2)[0]->getMatrixIndex()) + initTransM;
	}

	for(unsigned int i=1; i< s1->size(); i++){
		avp += exp(previous->getValueAt(k,l));
		if((*s2)[i]->isIsGap()){
			//Insert
			k++;
			if (previous == X)
				lnl += X->getTransitionProbabilityFromInsert();
			else if (previous == Y)
				lnl += X->getTransitionProbabilityFromDelete();
			else
				lnl+=X->getTransitionProbabilityFromMatch();
			lnl += ptmatrix->getLogEquilibriumFreq((*s1)[i]->getMatrixIndex());
			if (post) DUMP("I " << i << "\tlnl\t" << lnl << "\tposterior\t" << exp(previous->getValueAt(k,l)));
			previous = X;
		}
		else if((*s1)[i]->isIsGap()){
			//Delete
			l++;
			if (previous == X)
				lnl += Y->getTransitionProbabilityFromInsert();
			else if (previous == Y)
				lnl += Y->getTransitionProbabilityFromDelete();
			else
				lnl+=Y->getTransitionProbabilityFromMatch();
			lnl += ptmatrix->getLogEquilibriumFreq((*s2)[i]->getMatrixIndex());
			if (post) DUMP("D " << i << "\tlnl\t" << lnl << "\tposterior\t" << exp(previous->getValueAt(k,l)));
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
			lnl += ptmatrix->getLogPairTransition((*s1)[i]->getMatrixIndex(), (*s2)[i]->getMatrixIndex());
			posteriors[k][l] += 1;
			if (post) DUMP("M " << i << "\tlnl\t" << lnl << "\tposterior\t" << exp(previous->getValueAt(k,l)));
			previous = M;
		}

	}
	if (post) DUMP("Average posterior prob: " << (avp/s1->size()));
	//cerr << endl;
	return lnl;

}

double BackwardPairHMM::runAlgorithm()
{
	int i;
	int j;

	//initial
	double mL = Definitions::minMatrixLikelihood;
	//initial insert likelihood
	double xL = Definitions::minMatrixLikelihood;
	//initial delete likelihood
	double yL = Definitions::minMatrixLikelihood;

	double sX,sY,sM, sS;

	double bm, bx, by;

	double bmp, bxp,byp;

	M->initializeData(true);
	X->initializeData(true);
	Y->initializeData(true);

	for (j = ySize-1, i=xSize-1; j >= 0; j--)
	{
		bxp = (i==xSize-1) ? xL : X->getValueAt(i+1,j) + ptmatrix->getLogEquilibriumFreq((*seq1)[i]->getMatrixIndex());
		byp = (j==ySize-1) ? yL : Y->getValueAt(i,j+1) + ptmatrix->getLogEquilibriumFreq((*seq2)[j]->getMatrixIndex());
		bmp = (i==xSize-1 ||j==ySize-1) ? mL : M->getValueAt(i+1,j+1) + ptmatrix->getLogPairTransition((*seq1)[i]->getMatrixIndex(), (*seq2)[j]->getMatrixIndex());

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

			X->setValueAt(i, j, this->piI);
			Y->setValueAt(i, j, this->piD);
			M->setValueAt(i, j, this->piM);

/*
			X->setValueAt(i, j, 0);
			Y->setValueAt(i, j, 0);
			M->setValueAt(i, j, 0);
*/
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
		bxp = (i==xSize-1) ? xL : X->getValueAt(i+1,j) + ptmatrix->getLogEquilibriumFreq((*seq1)[i]->getMatrixIndex());
		byp = (j==ySize-1) ? yL : Y->getValueAt(i,j+1) + ptmatrix->getLogEquilibriumFreq((*seq2)[j]->getMatrixIndex());
		bmp = (i==xSize-1 ||j==ySize-1) ? mL : M->getValueAt(i+1,j+1) + ptmatrix->getLogPairTransition((*seq1)[i]->getMatrixIndex(), (*seq2)[j]->getMatrixIndex());

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


			X->setValueAt(i, j, this->piI);
			Y->setValueAt(i, j, this->piD);
			M->setValueAt(i, j, this->piM);

/*
			X->setValueAt(i, j, 0);
			Y->setValueAt(i, j, 0);
			M->setValueAt(i, j, 0);
*/
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

			bxp = X->getValueAt(i+1,j) + ptmatrix->getLogEquilibriumFreq((*seq1)[i]->getMatrixIndex());
			byp = Y->getValueAt(i,j+1) + ptmatrix->getLogEquilibriumFreq((*seq2)[j]->getMatrixIndex());
			bmp = M->getValueAt(i+1,j+1) + ptmatrix->getLogPairTransition((*seq1)[i]->getMatrixIndex(), (*seq2)[j]->getMatrixIndex());

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

	//FIXME - set to zero or ml, xl etc..
	sM = M->getValueAt(0, 0);
	sX = X->getValueAt(0, 0);
	sY = Y->getValueAt(0, 0);
	sS = maths->logSum(sM,sX,sY);

	DUMP("Backward results:");
	DUMP(" sX, sY, sM, sS " << sX << "\t" << sY << "\t" << sM << "\t" << sS);

	return sS* -1.0;
}

void BackwardPairHMM::calculateMaximumPosteriorMatrix() {
	//any state type will do
	this->MPstate = new PairwiseHmmMatchState(xSize,ySize);

	double tmpMax;
	unsigned int k,l;
	//no need to initialize data
	//set the P00 to 1 (ln(1) = 0)
	MPstate->setValueAt(0,0,0);
	for (unsigned int i = 1 ; i < xSize; i++)
		for (unsigned int j = 1; j < ySize; j++){
			k = i-1;
			l = j-1;
			tmpMax = std::max(MPstate->getValueAt(k,l)+M->getValueAt(i,j), std::max(MPstate->getValueAt(k,j)+X->getValueAt(i,j), MPstate->getValueAt(i,l)+Y->getValueAt(i,j)));
			MPstate->setValueAt(i,j,tmpMax);
		}

	DUMP("MPSTATE matrix ");
	dynamic_cast<DpMatrixFull*>(MPstate->getDpMatrix())->outputValues(0);


}

pair<string, string> BackwardPairHMM::getMPAlignment() {
	//tarceback through MPstate matrix
	//similar to viterbi traceback !
	DUMP("Backward HMM get MP alignment");
	pair<string, string> alignment;

	//unsigned char gapElem = this->substModel->getMatrixSize();

	//reserve memory for out strings (20% of gaps should be ok)
	alignment.first.reserve(max(xSize,ySize)*1.2);
	alignment.second.reserve(max(xSize,ySize)*1.2);


	double tm, ti, td;

	unsigned int i = xSize-1;
	unsigned int j = ySize-1;

	while(i > 0 && j > 0)
	{
		tm = MPstate->getValueAt(i-1,j-1);
		ti = MPstate->getValueAt(i-1,j);
		td = MPstate->getValueAt(i,j-1);

		if (tm >= ti && td <= ti){
			//MATCH
			alignment.first += (*seq1)[i-1]->getSymbol();
			alignment.second += (*seq2)[j-1]->getSymbol();
			i--;
			j--;
		}
		else if (ti >= td){
			//INSERT
			alignment.first += (*seq1)[i-1]->getSymbol();
			alignment.second += '-';//gapElem;
			i--;
		}
		else{
			//DELETE
			alignment.first += '-';//gapElem;
			alignment.second += (*seq2)[j-1]->getSymbol();
			j--;
		}
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


} /* namespace EBC */


