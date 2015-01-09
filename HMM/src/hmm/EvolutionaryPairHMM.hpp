/*
 * EvolutionaryPairHMM.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#ifndef EVOLUTIONARYPAIRHMM_HPP_
#define EVOLUTIONARYPAIRHMM_HPP_


#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/PMatrixDouble.hpp"
#include "core/TransitionProbabilities.hpp"
#include "core/Sequences.hpp"
#include "core/Maths.hpp"
#include "core/SequenceElement.hpp"

#include "hmm/PairwiseHmmStateBase.hpp"
#include "hmm/PairwiseHmmInsertState.hpp"
#include "hmm/PairwiseHmmDeleteState.hpp"
#include "hmm/PairwiseHmmMatchState.hpp"
#include "hmm/DpMatrixLoMem.hpp"

#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/IndelModel.hpp"
#include "models/SubstitutionModelBase.hpp"

#include "heuristics/Band.hpp"

namespace EBC
{

//Thats a 3 state pair-HMM with an insertion,deletion and match state
//The states are connected by the silent (blank) state
class EvolutionaryPairHMM
{
protected:

	Dictionary* dict;
	SubstitutionModelBase* substModel;

	PMatrixDouble* ptmatrix;

	IndelModel* indelModel;
	Sequences* inputSequences;
	TransitionProbabilities* tpb;
	Maths* maths;

	double * mlParameters;
	unsigned int totalParameters;
	unsigned int substParameters;
	unsigned int indelParameters;

	unsigned int xSize, ySize;

	//Bound scale
	unsigned int bandFactor;
	unsigned int bandSpan;
	unsigned int gammaRateCategories;

	bool bandingEnabled;

	bool equilibriumFreqs;

	Band* band;

	double initTransM;
	double initTransX;
	double initTransY;

	vector<SequenceElement> seq1;
	vector<SequenceElement> seq2;
	vector<SequenceElement>::iterator itS1, itS2;
	
	//cumulative likelihood for all 3 matrices
	double totalLikelihood;

	//state transition matrix
	double md[Definitions::stateCount][Definitions::stateCount];

	//gap probs;
	double e,g;

	//state equilibruim frequencies
	double piM, piI, piD;

	//the following assumes a fix HMM structure
	virtual void setTransitionProbabilities();

	virtual void calculateModels();

	virtual void initializeStates(Definitions::DpMatrixType mt);

	void getStateEquilibriums();


public:

	//Match state
	PairwiseHmmStateBase* M;
	//Insert state
	PairwiseHmmStateBase* X;
	//Delete state
	PairwiseHmmStateBase* Y;

	EvolutionaryPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, SubstitutionModelBase* smdl,
			IndelModel* imdl, Definitions::DpMatrixType, Band* bandObj, bool useEquilibriumFreqs);

	virtual ~EvolutionaryPairHMM();

	virtual double runAlgorithm()=0;

	void summarize();

	void setDivergenceTime(double time);

	unsigned int getIndelParameterCount()
	{
		return this->indelModel->getParamsNumber();
	}

	unsigned int getSubstitutionParameterCount()
	{
		return this->substModel->getParamsNumber();
	}

	unsigned int getTotalParameters() const
	{
		return totalParameters;
	}

	PairwiseHmmStateBase* getM()
	{
		return M;
	}

	PairwiseHmmStateBase* getX()
	{
		return X;
	}

	PairwiseHmmStateBase* getY()
	{
		return Y;
	}

	double getAlignmentLikelihood(vector<SequenceElement> s1, vector<SequenceElement> s2);

	double getTotalLikelihood() const {
		return totalLikelihood;
	}

	void setTotalLikelihood(double totalLikelihood) {
		this->totalLikelihood = totalLikelihood;
	}
};

} /* namespace EBC */
#endif /* EVOLUTIONARYPAIRHMM_HPP_ */
