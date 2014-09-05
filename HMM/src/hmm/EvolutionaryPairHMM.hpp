/*
 * EvolutionaryPairHMM.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#ifndef EVOLUTIONARYPAIRHMM_HPP_
#define EVOLUTIONARYPAIRHMM_HPP_

#include "models/SubstitutionModelBase.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "models/IndelModel.hpp"
#include "core/Sequences.hpp"
#include "core/Maths.hpp"
#include "hmm/PairwiseHmmStateBase.hpp"
#include "core/SequenceElement.hpp"
#include "hmm/PairwiseHmmInsertState.hpp"
#include "hmm/PairwiseHmmDeleteState.hpp"
#include "hmm/PairwiseHmmMatchState.hpp"
#include "hmm/DpMatrixLoMem.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"

namespace EBC
{

//Thats a 3 state pair-HMM with an insertion,deletion and match state
//The states are connected by the silent (blank) state
class EvolutionaryPairHMM
{
protected:

	Dictionary* dict;
	SubstitutionModelBase* substModel;
	IndelModel* indelModel;
	Sequences* inputSequences;
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

	vector<SequenceElement> seq1;
	vector<SequenceElement> seq2;
	vector<SequenceElement>::iterator itS1, itS2;

	//Match state
	PairwiseHmmStateBase* M;
	//Insert state
	PairwiseHmmStateBase* X;
	//Delete state
	PairwiseHmmStateBase* Y;

	inline void getBandWidth()
	{
		this->bandSpan = ySize/(bandFactor);
		DEBUG("Band span " << bandSpan);
	}

	inline bool withinBand(unsigned int line, int position, unsigned int width)
	{
		int low = line - width;
		int high = line + width;
		bool result = ((position >= low) && (position <= high));
		return result;
	}


	//the following assumes a fix HMM structure
	virtual void setTransitionProbabilities();

	virtual void calculateModels();

	virtual void initializeStates(Definitions::DpMatrixType mt);

	double* generateInitialSubstitutionParameters();

	double* generateInitialIndelParameters();

	double generateInitialDistanceParameter();


public:

	EvolutionaryPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dict, unsigned int rateCategories, Maths*,
			Definitions::ModelType model, bool banding, unsigned int bandPercentage, Definitions::DpMatrixType mt);

	virtual ~EvolutionaryPairHMM();

	virtual double runAlgorithm()=0;

	void summarize();

	void setModelParameters(std::vector<double> indel_params, std::vector<double> subst_params,double evolDistance, double alpha);

	void setModelFrequencies(double*);

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
};

} /* namespace EBC */
#endif /* EVOLUTIONARYPAIRHMM_HPP_ */
