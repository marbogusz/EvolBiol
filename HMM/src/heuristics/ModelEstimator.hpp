/*
 * ModelEstimator.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef MODELESTIMATOR_HPP_
#define MODELESTIMATOR_HPP_


#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include "core/HmmException.hpp"
#include "core/PMatrixTriple.hpp"

#include "models/SubstitutionModelBase.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/GTRModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"

#include "heuristics/ModelEstimator.hpp"
#include "heuristics/GuideTree.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletAligner.hpp"
#include "heuristics/StateTransitionEstimator.hpp"
#include "heuristics/SubstitutionModelEstimator.hpp"

#include "hmm/ViterbiPairHMM.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "hmm/BackwardPairHMM.hpp"
#include "hmm/DpMatrixFull.hpp"


#include <sstream>
#include <vector>
#include <array>

using namespace std;


namespace EBC
{

class ModelEstimator
{
protected:


	Dictionary* dict;

	SubstitutionModelBase* substModel;
	IndelModel* indelModel;

	Sequences* inputSequences;
	Maths* maths;

	GuideTree* gtree;
	TripletSamplingTree tst;

	StateTransitionEstimator* ste;
	SubstitutionModelEstimator* sme;

	TripletAligner* tal;
	ViterbiPairHMM* vphmm;

	double bestFwdAlpha;
	double bestFwdTm;

	bool estAlpha;
	bool estIndel;
	bool estSubst;

	vector<array<vector<unsigned char>*, 3> > tripleAlignments;
	vector<array<vector<unsigned char>*, 4> > pairAlignments;
	vector<array<vector<double>*, 2> > pairwisePosteriors;
	vector<array<unsigned int, 3> > tripletIdxs;
	vector<array<double, 3> > tripletDistances;
	vector<array<ForwardPairHMM*,2> > fwdHMMs;

	vector<double> substitutionParameters;
	vector<double> indelParameters;
	double alpha;

	unsigned int gammaRateCategories;

	unsigned int tripletIdxsSize;

	double userAlpha;

	bool estimateAlpha;

	void calculateInitialHMMs(Definitions::ModelType model);

	void estimateParameters();

	Definitions::ModelType model;

	vector<double> getInitialModelParameters();

public:
	ModelEstimator(Sequences* inputSeqs, Definitions::ModelType model,
			Definitions::OptimizationType ot,
			unsigned int rateCategories, double alpha, bool estimateAlpha);

	virtual ~ModelEstimator();



	vector<double> getSubstitutionParameters();
	vector<double> getIndelParameters();
	double getAlpha();

	void recalculateHMMs();

	GuideTree* getGuideTree()
	{
		return gtree;
	}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
