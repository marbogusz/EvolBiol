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

void BackwardPairHMM::getAlignment(HMMPathSample& sample)
{

	unsigned char gapElem = this->substModel->getMatrixSize(); // last matrix element is the gap ID!

	unsigned int i = xSize-1;
	unsigned int j = ySize-1;

	//choose initial state
	PairwiseHmmStateBase* currentState = NULL;
	PairwiseHmmStateBase* previousState = NULL;


	double tm, ti, td;

	tm = MPstate->getValueAt(i-1,j-1);
	ti = MPstate->getValueAt(i-1,j);
	td = MPstate->getValueAt(i,j-1);

	if (tm >= ti && tm >= td){
		previousState = M;
		sample.addSitePattern((*seq1)[i-1]->getMatrixIndex(),(*seq2)[j-1]->getMatrixIndex());
		i--;
		j--;
	}
	else if (ti >= td){
		//INSERT
		//DUMP("I");
		previousState = X;
		sample.addSitePattern((*seq1)[i-1]->getMatrixIndex(),gapElem);
		i--;
	}
	else{
		//DELETE
		//DUMP("D");
		previousState = Y;
		sample.addSitePattern(gapElem,(*seq2)[j-1]->getMatrixIndex());
		j--;
	}

	while(i > 0 && j > 0)
	{
		tm = MPstate->getValueAt(i-1,j-1);
		ti = MPstate->getValueAt(i-1,j);
		td = MPstate->getValueAt(i,j-1);

		//DUMP(tm << "\t" << ti << "\t" << td);

		if (tm >= ti && tm >= td){
			currentState = M;
			sample.addSitePattern((*seq1)[i-1]->getMatrixIndex(),(*seq2)[j-1]->getMatrixIndex());
			i--;
			j--;
		}
		else if (ti >= td){
			//INSERT
			//DUMP("I");
			currentState = X;
			sample.addSitePattern((*seq1)[i-1]->getMatrixIndex(),gapElem);
			i--;
		}
		else{
			//DELETE
			//DUMP("D");
			currentState = Y;
			sample.addSitePattern(gapElem,(*seq2)[j-1]->getMatrixIndex());
			j--;
		}
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
	DUMP("FORWARD MATRICES");
	dynamic_cast<DpMatrixFull*>(fwd->M->getDpMatrix())->outputValues(0);
	DUMP("BACKWARD MATRICES");
	dynamic_cast<DpMatrixFull*>(M->getDpMatrix())->outputValues(0);
	DUMP("\nINSERT");
	DUMP("FORWARD MATRICES");
	dynamic_cast<DpMatrixFull*>(fwd->X->getDpMatrix())->outputValues(0);
	DUMP("BACKWARD MATRICES");
	dynamic_cast<DpMatrixFull*>(X->getDpMatrix())->outputValues(0);

	DUMP("\nDELETE");
	DUMP("FORWARD MATRICES");
	dynamic_cast<DpMatrixFull*>(fwd->Y->getDpMatrix())->outputValues(0);
	DUMP("BACKWARD MATRICES");
	dynamic_cast<DpMatrixFull*>(Y->getDpMatrix())->outputValues(0);
	//DUMP("#####Match posteriors########");

	//DUMP("#####Insert posteriors########");

	//dynamic_cast<DpMatrixFull*>(X->getDpMatrix())->outputValues(0);
	//DUMP("#####Delete posteriors########");

	//dynamic_cast<DpMatrixFull*>(Y->getDpMatrix())->outputValues(0);

	DUMP("POSTERIORS MATRICES");
*/
	for (i = 1; i<=xSize-1; i++)
	{
		for (j = 1; j<=ySize-1; j++)
		{
			xval = X->getValueAt(i,j) + fwd->X->getValueAt(i,j) - fwdT;
			yval = Y->getValueAt(i,j) + fwd->Y->getValueAt(i,j) - fwdT;
			mval = M->getValueAt(i,j) + fwd->M->getValueAt(i,j) - fwdT;

			X->setValueAt(i,j,xval);
			Y->setValueAt(i,j,yval);
			M->setValueAt(i,j,mval);
		}
	}
/*
	DUMP("#####Match posteriors########");
	dynamic_cast<DpMatrixFull*>(M->getDpMatrix())->outputValues(0);
	DUMP("#####Insert posteriors########");
	dynamic_cast<DpMatrixFull*>(X->getDpMatrix())->outputValues(0);
	DUMP("#####Delete posteriors########");
	dynamic_cast<DpMatrixFull*>(Y->getDpMatrix())->outputValues(0);
*/
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

	double initProb = log(xi);

	M->initializeData(true);
	X->initializeData(true);
	Y->initializeData(true);

	//Last ROW
	for (j = ySize-1, i=xSize-1; j > 0; j--)
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
			X->setValueAt(i, j, initProb);
			Y->setValueAt(i, j, initProb);
			M->setValueAt(i, j, initProb);

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
	//LAST COLUMN
	for (i = xSize-1, j=ySize-1; i > 0; i--)
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

			X->setValueAt(i, j, initProb);
			Y->setValueAt(i, j, initProb);
			M->setValueAt(i, j, initProb);

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

	//FIRST INSERTION boundary
	X->setValueAt(xSize-1,0,ptmatrix->getLogEquilibriumFreq((*seq2)[0]->getMatrixIndex())+
			Y->getTransitionProbabilityFromInsert()+Y->getValueAt(xSize-1,1));

	//FIRST DELETION boundary
	Y->setValueAt(0,ySize-1,ptmatrix->getLogEquilibriumFreq((*seq1)[0]->getMatrixIndex())+
				X->getTransitionProbabilityFromDelete()+X->getValueAt(1,ySize-1));


	//EVERTYHING ELSE excl first X col and 1st Y row
	if(this->band == NULL){
		for (i = xSize-2; i > 0; i--)
		{
			for (j = ySize-2; j > 0; j--)
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

		//set final inserts and deletes!

		//first X col
		for (j=0,i=xSize-2; i > 0; i--){

			bxp = X->getValueAt(i+1,j) + ptmatrix->getLogEquilibriumFreq((*seq1)[i]->getMatrixIndex());
			byp = Y->getValueAt(i,j+1) + ptmatrix->getLogEquilibriumFreq((*seq2)[j]->getMatrixIndex());
			bmp = M->getValueAt(i+1,j+1) + ptmatrix->getLogPairTransition((*seq1)[i]->getMatrixIndex(), (*seq2)[j]->getMatrixIndex());

			bx = maths->logSum(M->getTransitionProbabilityFromInsert() + bmp,
					X->getTransitionProbabilityFromInsert() + bxp,
					Y->getTransitionProbabilityFromInsert() + byp);
			X->setValueAt(i, j, bx);
		}
		//first Y row
		for (i=0,j=ySize-2; j > 0; j--){
			bxp = X->getValueAt(i+1,j) + ptmatrix->getLogEquilibriumFreq((*seq1)[i]->getMatrixIndex());
			byp = Y->getValueAt(i,j+1) + ptmatrix->getLogEquilibriumFreq((*seq2)[j]->getMatrixIndex());
			bmp = M->getValueAt(i+1,j+1) + ptmatrix->getLogPairTransition((*seq1)[i]->getMatrixIndex(), (*seq2)[j]->getMatrixIndex());

			by = maths->logSum(M->getTransitionProbabilityFromDelete() + bmp,
								X->getTransitionProbabilityFromDelete() + bxp,
								Y->getTransitionProbabilityFromDelete() + byp);

			Y->setValueAt(i, j, by);
		}
	}
	else{
		int loD, hiD;
		int iMin = xSize-2;
		int iMax = 0;

		//FIXME - I make a simplifying assumption that the band is the same for all 3 matrices!
		//Use bracket D bracket to calculate, but the first row for M and I needs to be zeroed!

		for (j = ySize-2; j >= 0; j--){
			auto bracketD = band->getDeleteRangeAt(j);

			loD = max((bracketD.first), iMax);
			hiD = min(bracketD.second, iMin);


			for (int i = hiD; i >= loD; i--){
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
		//do the first column
		auto bracketI = band->getInsertRangeAt(0);
		int loI = max((bracketI.first), iMax);
		int hiI = min(bracketI.second, iMin);


		for (j=0,i=hiI; i >= loI; i--){
			bxp = X->getValueAt(i+1,j) + ptmatrix->getLogEquilibriumFreq((*seq1)[i]->getMatrixIndex());
			byp = Y->getValueAt(i,j+1) + ptmatrix->getLogEquilibriumFreq((*seq2)[j]->getMatrixIndex());
			bmp = M->getValueAt(i+1,j+1) + ptmatrix->getLogPairTransition((*seq1)[i]->getMatrixIndex(), (*seq2)[j]->getMatrixIndex());

			bx = maths->logSum(M->getTransitionProbabilityFromInsert() + bmp,
					X->getTransitionProbabilityFromInsert() + bxp,
					Y->getTransitionProbabilityFromInsert() + byp);
			X->setValueAt(i, j, bx);
		}

		//zero the first col
		M->getDpMatrix()->setWholeRow(0, Definitions::minMatrixLikelihood);
		X->getDpMatrix()->setWholeRow(0, Definitions::minMatrixLikelihood);
	}

	//FIXME - remove to improve performance by 0.000000001% :)
	bmp  = M->getValueAt(1,1) + ptmatrix->getLogPairTransition((*seq1)[0]->getMatrixIndex(), (*seq2)[0]->getMatrixIndex());
	bxp  = X->getValueAt(1,0) + ptmatrix->getLogEquilibriumFreq((*seq1)[0]->getMatrixIndex());
	byp  = Y->getValueAt(0,1) + ptmatrix->getLogEquilibriumFreq((*seq2)[0]->getMatrixIndex());
	bm = bmp + initTransM;
	bx = bxp + initTransX;
	by = byp + initTransY;
	M->setValueAt(0, 0, maths->logSum(bm,bx,by));
	sS = maths->logSum(bm,bx,by);

	//DUMP("Backward results:");
	//DUMP(" sX, sY, sM, sS " << sX << "\t" << sY << "\t" << sM << "\t" << sS);

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

	//TODO - we can band it as well

	for (unsigned int i = 1 ; i < xSize; i++)
		for (unsigned int j = 1; j < ySize; j++){
			k = i-1;
			l = j-1;
			tmpMax = std::max(MPstate->getValueAt(k,l)+M->getValueAt(i,j), std::max(MPstate->getValueAt(k,j)+X->getValueAt(i,j), MPstate->getValueAt(i,l)+Y->getValueAt(i,j)));
			MPstate->setValueAt(i,j,tmpMax);
		}

	//DUMP("MPSTATE matrix ");
	//dynamic_cast<DpMatrixFull*>(MPstate->getDpMatrix())->outputValues(0);


}
pair<vector<double>*, pair<vector<unsigned char>*, vector<unsigned char>*> >
BackwardPairHMM::getMPDWithPosteriors(){
	DUMP("Backward HMM get MPD alignment with posteriors");
	pair<vector<double>*, pair<vector<unsigned char>*, vector<unsigned char>*> >
	ret = make_pair(new vector<double>(),make_pair(new vector<unsigned char>(), new vector<unsigned char>()));

	unsigned char gapElem = this->substModel->getMatrixSize(); // last matrix element is the gap ID!

	double tm, ti, td;

	unsigned int i = xSize-1;
	unsigned int j = ySize-1;

	//FIXME - which index should I consider while performing a readback ???

	while(i > 0 && j > 0)
	{
		tm = MPstate->getValueAt(i-1,j-1);
		ti = MPstate->getValueAt(i-1,j);
		td = MPstate->getValueAt(i,j-1);

		//DUMP(tm << "\t" << ti << "\t" << td);

		if (tm >= ti && tm >= td){
			//MATCH
			//DUMP("M");
			ret.second.first->push_back((*seq1)[i-1]->getMatrixIndex());
			ret.second.second->push_back((*seq2)[j-1]->getMatrixIndex());
			ret.first->push_back(M->getValueAt(i,j));
			i--;
			j--;

		}
		else if (ti >= td){
			//INSERT
			//DUMP("I");
			ret.second.first->push_back((*seq1)[i-1]->getMatrixIndex());
			ret.second.second->push_back(gapElem);//gapElem;
			ret.first->push_back(X->getValueAt(i,j));
			i--;

		}
		else{
			//DELETE
			//DUMP("D");
			ret.second.first->push_back(gapElem);//gapElem;
			ret.second.second->push_back((*seq2)[j-1]->getMatrixIndex());
			ret.first->push_back(Y->getValueAt(i,j));
			j--;

		}
	}


	if (j==0)
	{
		while(i > 0){
			ret.second.first->push_back((*seq1)[i-1]->getMatrixIndex());
			ret.second.second->push_back(gapElem);
			ret.first->push_back(X->getValueAt(i,j));
			i--;

		}
	}
	else if (i==0)
	{
		while(j > 0){
			ret.second.second->push_back((*seq2)[j-1]->getMatrixIndex());
			ret.second.first->push_back(gapElem);
			ret.first->push_back(Y->getValueAt(i,j));
			j--;
		}
	}
	//deal with the last row or column

	reverse(ret.second.first->begin(), ret.second.first->end());
	reverse(ret.second.second->begin(), ret.second.second->end());
	reverse(ret.first->begin(), ret.first->end());
	return ret;
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

		//DUMP(tm << "\t" << ti << "\t" << td);

		if (tm >= ti && tm >= td){
			//MATCH
			//DUMP("M");
			alignment.first += (*seq1)[i-1]->getSymbol();
			alignment.second += (*seq2)[j-1]->getSymbol();
			i--;
			j--;
		}
		else if (ti >= td){
			//INSERT
			//DUMP("I");
			alignment.first += (*seq1)[i-1]->getSymbol();
			alignment.second += '-';//gapElem;
			i--;
		}
		else{
			//DELETE
			//DUMP("D");
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


