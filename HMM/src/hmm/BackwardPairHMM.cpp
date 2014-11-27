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
#include <random>


namespace EBC
{


BackwardPairHMM::BackwardPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, SubstitutionModelBase* smdl,
		IndelModel* imdl,  Definitions::DpMatrixType mt ,Band* bandObj) : dpMatrixCalculated(false), posteriorsCalculated(false),
		EvolutionaryPairHMM(s1,s2, smdl, imdl, mt, bandObj)
{

}

BackwardPairHMM::~BackwardPairHMM()
{
}

void BackwardPairHMM::calculatePosteriors(ForwardPairHMM* fwd)
{
	if (!dpMatrixCalculated)
		throw HmmException("Error - attempting to calculate posterior probabilities without running backward algorithm first");

	DEBUG("Calculating posterior probabilities");

	int i,j;
	double xval, yval, mval;
	//totalForward
	double fwdT;

	fwdT = fwd->getTotalLikelihood();

	for (i = 0; i<xSize-1; i++)
	{
		for (j = 0; j<ySize-1; j++)
		{
			xval = X->getValueAt(i,j) + fwd->X->getValueAt(i,j) - fwdT;
			if (xval > 0.01){
				ERROR("Posterior lnl > 0 X " << xval << " bw " << X->getValueAt(i,j) << " fwx " <<  fwd->X->getValueAt(i,j) << " fwtM "
						<<  fwdT << " fwtX" << fwd->X->getValueAt(xSize-1,ySize-1) << " fwtY" << fwd->Y->getValueAt(xSize-1,ySize-1));
			}
			yval = Y->getValueAt(i,j) + fwd->Y->getValueAt(i,j) - fwdT;
			if (yval > 0.01){
				ERROR("Posterior lnl > 0 Y " << yval << " bw " << Y->getValueAt(i,j) << " fwy " <<  fwd->Y->getValueAt(i,j) << " fwtM "
						<<  fwdT << " fwtX" << fwd->X->getValueAt(xSize-1,ySize-1) << " fwtY" << fwd->Y->getValueAt(xSize-1,ySize-1));
			}
			mval = M->getValueAt(i,j) + fwd->M->getValueAt(i,j) - fwdT;
			if (mval > 0.01){
				ERROR("Posterior lnl > 0 M " << mval << " bw " << M->getValueAt(i,j) << " fwm " <<  fwd->M->getValueAt(i,j) << " fwtM "
						<<  fwdT << " fwtX" << fwd->X->getValueAt(xSize-1,ySize-1) << " fwtY" << fwd->Y->getValueAt(xSize-1,ySize-1));
			}

			X->setValueAt(i,j,xval);
			Y->setValueAt(i,j,yval);
			M->setValueAt(i,j,mval);
		}
	}

	posteriorsCalculated = true;

	DUMP("Match");
	dynamic_cast<DpMatrixFull*>(M->getDpMatrix())->outputValues(0);
	DUMP("Insert");
	dynamic_cast<DpMatrixFull*>(X->getDpMatrix())->outputValues(0);
	DUMP("Delete");
	dynamic_cast<DpMatrixFull*>(Y->getDpMatrix())->outputValues(0);

/*
	for (int i=0;i<10;i++)
		for(int j=0;j<10;j++)
		{
			double mt,in,dl;
			//DUMP("M-H " << exp(M->getValueAt(i,j+1)) << " M-V " << exp(M->getValueAt(i+1,j)) << " M-D " << exp(M->getValueAt(i+1,j+1)));
			//DUMP("I-H " << exp(X->getValueAt(i,j+1)) << " I-V " << exp(X->getValueAt(i+1,j)) << " I-D " << exp(X->getValueAt(i+1,j+1)));
			//DUMP("D-H " << exp(Y->getValueAt(i,j+1)) << " D-V " << exp(Y->getValueAt(i+1,j)) << " D-D " << exp(Y->getValueAt(i+1,j+1)));
			//DUMP("M-D + I-V + D-H " <<  (exp(M->getValueAt(i+1,j+1))+exp(X->getValueAt(i+1,j))+exp(Y->getValueAt(i,j+1))));
			mt = exp(M->getValueAt(i,j+1)) + exp(M->getValueAt(i+1,j)) + exp(M->getValueAt(i+1,j+1));
			in = exp(X->getValueAt(i,j+1)) + exp(X->getValueAt(i+1,j)) + exp(X->getValueAt(i+1,j+1));
			dl = exp(Y->getValueAt(i,j+1)) + exp(Y->getValueAt(i+1,j)) + exp(Y->getValueAt(i+1,j+1));
			DUMP(" M " << mt << " I " << in << " D " << dl << "\t\t\t Total " << (mt+in+dl) );
		}
*/
}

pair<string, string>& BackwardPairHMM::sampleAlignment(string&seq_a, string& seq_b)
{
	//FIXME - reference from stack
	if (!posteriorsCalculated)
		throw HmmException("Error - attempting to sample an alignment  without calculating posterior probabilities");

	cerr << seq_a << endl;
	cerr << seq_b << endl;

	cerr << xSize << endl;
	cerr << ySize << endl;

	string s1;
	string s2;
	s1.reserve(max(xSize,ySize)*1.2);
	s2.reserve(max(xSize,ySize)*1.2);
	pair<string, string> alignment = make_pair(s1,s2);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(1, 2);

	unsigned i = 0;
	unsigned j = 0;
	double mt,in,dl, rnbr;

	while(i < xSize-1 || j < ySize-1)
	{
		//get probabilities
		mt = exp(M->getValueAt(i,j+1)) + exp(M->getValueAt(i+1,j)) + exp(M->getValueAt(i+1,j+1));
		in = exp(X->getValueAt(i,j+1)) + exp(X->getValueAt(i+1,j)) + exp(X->getValueAt(i+1,j+1));
		dl = exp(Y->getValueAt(i,j+1)) + exp(Y->getValueAt(i+1,j)) + exp(Y->getValueAt(i+1,j+1));

		std::uniform_real_distribution<> dis(0, mt+in+dl);
		rnbr = dis(gen);
		cerr << " M " << mt << " I " << in << " D " << dl << "Sampled " << rnbr;
		if (rnbr <= mt){
			alignment.first += seq_a[i];
			alignment.second += seq_b[j];
			i++;
			j++;
		}

		else if(rnbr <= mt+in){
			alignment.first += seq_a[i];
			alignment.second += '-';
			i++;
		}
		else{
			alignment.first += '-';
			alignment.second += seq_b[j];
			j++;
		}
		cerr << alignment.first << i << endl;
		cerr << alignment.second << j << endl;
	}

	cerr << "DONE" << endl;
	return alignment;
}

double BackwardPairHMM::runAlgorithm()
{
	calculateModels();
	setTransitionProbabilities();

	int i;
	int j;

	double sX,sY,sM, sS;

	double bm, bx, by;

	double bmp, bxp,byp;

	M->initializeData(true);
	X->initializeData(true);
	Y->initializeData(true);


	for (j = ySize-1, i=xSize-1; j >= 0; j--)
	{
		bxp = (i==xSize-1) ? -100000 : X->getValueAt(i+1,j) + ptmatrix->getLogEquilibriumFreq(seq1[i].getMatrixIndex());
		byp = (j==ySize-1) ? -100000 : Y->getValueAt(i,j+1) + ptmatrix->getLogEquilibriumFreq(seq2[j].getMatrixIndex());
		bmp = (i==xSize-1 ||j==ySize-1) ? -100000 : M->getValueAt(i+1,j+1) + ptmatrix->getLogPairTransition(seq1[i].getMatrixIndex(), seq2[j].getMatrixIndex());

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
			X->setValueAt(i, j, 0);
			Y->setValueAt(i, j, 0);
			M->setValueAt(i, j, 0);
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
		bxp = (i==xSize-1) ? -100000 : X->getValueAt(i+1,j) + ptmatrix->getLogEquilibriumFreq(seq1[i].getMatrixIndex());
		byp = (j==ySize-1) ? -100000 : Y->getValueAt(i,j+1) + ptmatrix->getLogEquilibriumFreq(seq2[j].getMatrixIndex());
		bmp = (i==xSize-1 ||j==ySize-1) ? -100000 : M->getValueAt(i+1,j+1) + ptmatrix->getLogPairTransition(seq1[i].getMatrixIndex(), seq2[j].getMatrixIndex());

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
			X->setValueAt(i, j, 0);
			Y->setValueAt(i, j, 0);
			M->setValueAt(i, j, 0);
		}
		else
		{
			X->setValueAt(i, j, bx);
			Y->setValueAt(i, j, by);
			M->setValueAt(i, j, bm);
		}
	}

	if(this->band == NULL)
	{
		for (i = xSize-2; i >= 0; i--)
		{
			for (j = ySize-2; j >= 0; j--)
			{

				bxp = X->getValueAt(i+1,j) +   ptmatrix->getLogEquilibriumFreq(seq1[i].getMatrixIndex());
				byp = Y->getValueAt(i,j+1) +   ptmatrix->getLogEquilibriumFreq(seq2[j].getMatrixIndex());
				bmp = M->getValueAt(i+1,j+1) + ptmatrix->getLogPairTransition (seq1[i].getMatrixIndex(), seq2[j].getMatrixIndex());

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
	}
	else
	{
		int loI, hiI, loD, hiD, loM, hiM;
		//banding!
		for (j = ySize-2; j >= 0; j--)
		{
			//FIXME - range should be a reference perhaps
			auto bracketM = band->getMatchRangeAt(j);

			//Only using the M range as backwards banding should not be used for custom bands.
			//Only the default one

			hiM = bracketM.first;
			hiM = hiM > xSize-2 ? xSize-2:hiM;
			if (hiM != -1)
			{
				hiM = hiM > xSize-2 ? xSize-2:hiM;
				for(i=hiM; i>=loM; i--)
				{

					bxp = X->getValueAt(i+1,j) +   ptmatrix->getLogEquilibriumFreq(seq1[i].getMatrixIndex());
					byp = Y->getValueAt(i,j+1) +   ptmatrix->getLogEquilibriumFreq(seq2[j].getMatrixIndex());
					bmp = M->getValueAt(i+1,j+1) + ptmatrix->getLogPairTransition (seq1[i].getMatrixIndex(), seq2[j].getMatrixIndex());

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
		}
	}
	sM = M->getValueAt(0, 0);
	sX = X->getValueAt(0, 0);
	sY = Y->getValueAt(0, 0);
	sS = maths->logSum(sM,sX,sY);

	this->setTotalLikelihood(sS);

	DUMP("Backward results:");
	DUMP(" sX, sY, sM, sS " << sX << "\t" << sY << "\t" << sM << "\t" << sS);

	dpMatrixCalculated = true;

	return sS* -1.0;
}


} /* namespace EBC */


