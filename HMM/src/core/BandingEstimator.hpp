/*
 * BandingEstimator.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef BANDINGESTIMATOR_HPP_
#define BANDINGESTIMATOR_HPP_


#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/BandingEstimator.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include "core/HmmException.hpp"
#include "core/Optimizer.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "hmm/ViterbiPairHMM.hpp"
#include "models/SubstitutionModelBase.hpp"
#include "models/IndelModel.hpp"
#include "core/IOptimizable.hpp"
#include "heuristics/GuideTree.hpp"
#include "heuristics/BandCalculator.hpp"

#include <vector>
#include <sstream>

using namespace std;


namespace EBC
{

class BandingEstimator : public IOptimizable
{

private:

protected:

	Optimizer* bfgs;
	Dictionary* dict;
	SubstitutionModelBase* substModel;
	IndelModel* indelModel;
	Sequences* inputSequences;
	Maths* maths;
	GuideTree* gt;

	unsigned int bandFactor;
	unsigned int bandSpan;
	unsigned int gammaRateCategories;

	bool bandingEnabled;

	bool estimateSubstitutionParams;
	bool estimateIndelParams;
	bool estimateDivergence;
	bool estimateAlpha;

	unsigned int pairCount;

	vector<EvolutionaryPairHMM*> hmms;

	OptimizedModelParameters* modelParams;

public:
	BandingEstimator(Definitions::AlgorithmType at, Sequences* inputSeqs, Definitions::ModelType model,std::vector<double> indel_params,
			std::vector<double> subst_params, Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, GuideTree* gt);

	virtual ~BandingEstimator();

	void runIteration(const column_vector& bfgsParameters);

	double runIteration();

	void outputDistanceMatrix(stringstream&);

	vector<double> getOptimizedTimes()
	{
		return this->modelParams->getDivergenceTimes();
	}

	//ModelParameters getMlParameters()
	//{
	//	return this->modelParameters;
	//}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
