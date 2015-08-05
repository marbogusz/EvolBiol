/*
 * BandingEstimator.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef BANDINGESTIMATOR_HPP_
#define BANDINGESTIMATOR_HPP_



#include "core/IOptimizable.hpp"
#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/BandingEstimator.hpp"
#include "core/DistanceMatrix.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include "core/HmmException.hpp"
#include "core/Optimizer.hpp"
#include "core/BrentOptimizer.hpp"
#include "core/PairHmmCalculationWrapper.hpp"

#include "models/SubstitutionModelBase.hpp"
#include "models/IndelModel.hpp"

#include "heuristics/GuideTree.hpp"
#include "heuristics/BandCalculator.hpp"
#include "heuristics/Band.hpp"

#include "hmm/ForwardPairHMM.hpp"
#include "hmm/ViterbiPairHMM.hpp"

#include <vector>
#include <sstream>

using namespace std;


namespace EBC
{

class BandingEstimator : public IOptimizable
{

private:

protected:

	BrentOptimizer* numopt;
	Dictionary* dict;
	SubstitutionModelBase* substModel;
	IndelModel* indelModel;
	Sequences* inputSequences;
	Maths* maths;
	GuideTree* gt;

	Definitions::AlgorithmType algorithm;

	unsigned int bandFactor;
	unsigned int bandSpan;
	unsigned int gammaRateCategories;

	bool bandingEnabled;

	bool estimateSubstitutionParams;
	bool estimateIndelParams;
	bool estimateDivergence;
	bool estimateAlpha;

	unsigned int pairCount;

	//vector<EvolutionaryPairHMM*> hmms;
	//delete bands in the destructor
	//vector<Band*> bands;
	vector<double> divergenceTimes;

	OptimizedModelParameters* modelParams;

public:
	BandingEstimator(Definitions::AlgorithmType at, Sequences* inputSeqs, Definitions::ModelType model,std::vector<double> indel_params,
			std::vector<double> subst_params, Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, GuideTree* gt);

	virtual ~BandingEstimator();

	void runIteration(const column_vector& bfgsParameters);

	double runIteration();

	void outputDistanceMatrix(stringstream&);

	void optimizePairByPair();

	vector<double> getOptimizedTimes()
	{
		return this->divergenceTimes;
	}

	//ModelParameters getMlParameters()
	//{
	//	return this->modelParameters;
	//}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
