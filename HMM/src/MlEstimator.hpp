/*
 * MlEstimator.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef MLESTIMATOR_HPP_
#define MLESTIMATOR_HPP_


#include "OptimizedModelParameters.hpp"
#include "Definitions.hpp"
#include "MlEstimator.hpp"
#include "ForwardPairHMM.hpp"
#include "ViterbiPairHMM.hpp"
#include "SubstitutionModelBase.hpp"
#include "Maths.hpp"
#include "Dictionary.hpp"
#include "IndelModel.hpp"
#include "Sequences.hpp"
#include <sstream>
#include "HmmException.hpp"

#include <vector>

using namespace std;


namespace EBC
{

class MlEstimator
{

private:

	//BFGS optimization wrapper for dlib
	class BFGS
	{
	protected:
		column_vector initParams;
		column_vector lowerBounds;
		column_vector upperBounds;

		unsigned int paramsCount;

		MlEstimator* parent;

		Definitions::OptimizationType optimizationType;

	public:
		BFGS(MlEstimator* enclosing, Definitions::OptimizationType ot);
		virtual ~BFGS();
		void optimize();

		double objectiveFunction(const column_vector& m);

		const column_vector objectiveFunctionDerivative(const column_vector& m);
	};


protected:

	BFGS* bfgs;

	Dictionary* dict;
	SubstitutionModelBase* substModel;
	IndelModel* indelModel;
	Sequences* inputSequences;
	Maths* maths;

	unsigned int gammaRateCategories;

	bool estimateSubstitutionParams;
	bool estimateIndelParams;
	bool estimateDivergence;
	bool estimateAlpha;

	unsigned int pairCount;

	//for viterbi calculation
	vector<ViterbiPairHMM*> hmms;

	//
	vector<SubstitutionModelBase*> substs;

	OptimizedModelParameters* modelParams;

	bool useViterbi;

public:
	MlEstimator(Sequences* inputSeqs, Definitions::ModelType model,std::vector<double> indel_params,
			std::vector<double> subst_params,Definitions::OptimizationType ot,
			unsigned int rateCategories, double alpha, bool estimateAlpha, double userTime, bool useViterbi);

	virtual ~MlEstimator();

	void runIteration(const column_vector& bfgsParameters);

	double runIteration();

	vector<double>& getOptimizedTimes()
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
