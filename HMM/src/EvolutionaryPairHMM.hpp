/*
 * EvolutionaryPairHMM.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#ifndef EVOLUTIONARYPAIRHMM_HPP_
#define EVOLUTIONARYPAIRHMM_HPP_

#include "SubstitutionModel.hpp"
#include "Maths.hpp"
#include "Dictionary.hpp"
#include "IndelModel.hpp"
#include "Sequences.hpp"
#include "Maths.hpp"
#include "PairHmmStateBase.hpp"
#include "SequenceElement.hpp"

namespace EBC
{

//Thats a 3 state pair-HMM with an insertion,deletion and match state
//The states are connected by the silent (blank) state
class EvolutionaryPairHMM
{
protected:

	Dictionary* dict;
	SubstitutionModel* substModel;
	IndelModel* indelModel;
	Sequences* inputSequences;
	Maths* maths;

	double * mlParameters;
	unsigned int totalParameters;
	unsigned int substParameters;
	unsigned int indelParameters;

	unsigned int xSize, ySize;

	//TODO - remove after testing
	double testFreqs[4];

	vector<SequenceElement> seq1;
	vector<SequenceElement> seq2;
	vector<SequenceElement>::iterator itS1, itS2;

	//Match state
	PairHmmStateBase* M;
	//Insert state
	PairHmmStateBase* X;
	//Delete state
	PairHmmStateBase* Y;

	//the following assumes a fix HMM structure
	virtual void setTransitionProbabilities();

	virtual void initializeModels()=0;

	virtual void calculateModels();

	virtual void getSequencePair();

	virtual void initializeStates() = 0;



public:
	EvolutionaryPairHMM(Sequences* inputSeqs);

	virtual ~EvolutionaryPairHMM();

	void summarize();
};

} /* namespace EBC */
#endif /* EVOLUTIONARYPAIRHMM_HPP_ */
