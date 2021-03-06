/*
 * BandCalculator.hpp
 *
 *  Created on: Oct 27, 2014
 *      Author: root
 */

#ifndef HEURISTICS_BANDCALCULATOR_HPP_
#define HEURISTICS_BANDCALCULATOR_HPP_

#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/SubstitutionModelBase.hpp"

#include "core/PMatrixDouble.hpp"
#include "core/TransitionProbabilities.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Definitions.hpp"
#include "models/IndelModel.hpp"
#include "core/Sequences.hpp"

#include "hmm/ForwardPairHMM.hpp"
#include "hmm/BackwardPairHMM.hpp"

#include "heuristics/Band.hpp"

#include<vector>

namespace EBC
{

class BandCalculator
{
protected:

	vector<ForwardPairHMM*> fwd;
	BackwardPairHMM* bwd;

	vector<SequenceElement*>* seq1;
	vector<SequenceElement*>* seq2;

	SubstitutionModelBase* substModel;
	IndelModel* indelModel;

	double time;


	PMatrixDouble* ptMatrix;
	TransitionProbabilities* trProbs;

	Band* band;

	double posteriorLikelihoodLimit;
	double posteriorLikelihoodDelta;

	double bestTime;
	double accuracy;

	double leftBound;
	double rightBound;

	void processPosteriorProbabilities(BackwardPairHMM* hmm, Band* band);

public:
	BandCalculator(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, SubstitutionModelBase* sm, IndelModel* im, double divergenceTime);
	virtual ~BandCalculator();

	inline Band* getBand()
	{
		return this->band;
	}

	double getClosestDistance();

	double getBrentAccuracy();

	double getLeftBound() {
		return leftBound;
	}

	double getRightBound(){
		return rightBound;
	}
};

} /* namespace EBC */

#endif /* HEURISTICS_BANDCALCULATOR_HPP_ */
