/*
 * HMMEstimator.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef HMMESTIMATOR_HPP_
#define HMMESTIMATOR_HPP_


#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include "core/HmmException.hpp"
#include "core/PMatrixTriple.hpp"
#include "core/IOptimizable.hpp"
#include "core/Optimizer.hpp"

#include "models/SubstitutionModelBase.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/GTRModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"

#include "heuristics/GuideTree.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletAligner.hpp"

#include "hmm/ViterbiPairHMM.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "hmm/BackwardPairHMM.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "hmm/BackwardPairHMM.hpp"
#include "hmm/DpMatrixFull.hpp"

#include "sampling/PairHMMSampler.hpp"

#include <sstream>
#include <vector>
#include <array>

using namespace std;


namespace EBC
{

class HMMEstimator : public IOptimizable
{
protected:


	Dictionary* dict;

	SubstitutionModelBase* substModel;
	IndelModel* indelModel;

	Sequences* inputSequences;
	Maths* maths;

	GuideTree* gtree;
	TripletSamplingTree tst;

	TripletAligner* tal;
	ViterbiPairHMM* vphmm;

	OptimizedModelParameters* modelParams;
	Optimizer* bfgs;

	vector<array<unsigned int, 3> > tripletIdxs;

	vector<PairHMMSampler> sampleWorkers;

	unsigned int gammaRateCategories;

	unsigned int tripletIdxsSize;

	bool estimateAlpha;

	void calculateInitialPairs(Definitions::ModelType model);

	double runIteration();

	void optimise();

public:
	HMMEstimator(Sequences* inputSeqs, Definitions::ModelType model,
			Definitions::OptimizationType ot,
			unsigned int rateCategories, double alpha, bool estimateAlpha) ;

	virtual ~HMMEstimator();


	vector<double> getSubstitutionParameters();
	vector<double> getIndelParameters();

	double getAlpha();

	GuideTree* getGuideTree()
	{
		return gtree;
	}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
