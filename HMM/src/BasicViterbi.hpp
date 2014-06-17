/*
 * BasicViterbi.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#ifndef BASICVITERBI_HPP_
#define BASICVITERBI_HPP_

#include "SubstitutionModelBase.hpp"
#include "Maths.hpp"
#include "Dictionary.hpp"
#include "IndelModel.hpp"
#include "Sequences.hpp"
#include "Maths.hpp"
#include "PairwiseHmmStateBase.hpp"
#include "SequenceElement.hpp"

namespace EBC
{

//Thats a 3 state pair-HMM with an insertion,deletion and match state
//The states are connected by the silent (blank) state
class BasicViterbi
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

	vector<SequenceElement> seq1;
	vector<SequenceElement> seq2;
	vector<SequenceElement>::iterator itS1, itS2;

	//Match state
	PairwiseHmmStateBase* M;
	//Insert state
	PairwiseHmmStateBase* X;
	//Delete state
	PairwiseHmmStateBase* Y;

	//the following assumes a fix HMM structure
	void setTransitionProbabilities();

	void generateInitialParameters();

	void initializeModels();

	void calculateModels();

	void getSequencePair();

	void initializeStates();

	double getMax(double m, double x, double y, unsigned int i, unsigned int j, PairwiseHmmStateBase* state);

public:
	BasicViterbi(Sequences* inputSeqs, Definitions::ModelType model,std::vector<double> substParams, double distance, std::vector<double> indelParams,
			unsigned int rateCategories, double alpha, double* estimatedPArams = NULL);

	virtual ~BasicViterbi();

	void runViterbiAlgorithm();

	void getResults(stringstream&);
};

} /* namespace EBC */
#endif 
