/*
 * PairwiseEstimator.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef PAIRWISEESTIMATOR_HPP_
#define PAIRWISEESTIMATOR_HPP_


#include "OptimizedModelParameters.hpp"
#include <vector>

using namespace std;


namespace EBC
{

class PairwiseEstimator
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

		PairwiseEstimator* parent;

		Definitions::OptimizationType optimizationType;

	public:
		BFGS(PairwiseEstimator* enclosing, Definitions::OptimizationType ot);
		virtual ~BFGS();
		void optimize();

		double objectiveFunction(const column_vector& m);

		const column_vector objectiveFunctionDerivative(const column_vector& m);
	};



protected:


	BFGS* bfgs;

	unsigned int bandFactor;
	unsigned int bandSpan;
	unsigned int gammaRateCategories;

	double initialAlpha;

	bool bandingEnabled;

	bool estimateSubstitutionParams;
	bool estimateIndelParams;
	bool estimateDivergence;
	bool estimateAlpha;

	vector<double> userIndelParameters;
	vector<double> userSubstParameters;
	vector<double> divergenceTimes;

public:
	PairwiseEstimator(Sequences* inputSeqs, Definitions::ModelType model,std::vector<double> indel_params,
			std::vector<double> subst_params, Definitions::OptimizationType ot, bool banding,
			unsigned int bandPercentage, unsigned int rateCategories, double alpha, bool estimateAlpha);

	virtual ~PairwiseEstimator();

	void runIteration(const column_vector& bfgsParameters);

	double runOptimization();

	//ModelParameters getMlParameters()
	//{
	//	return this->modelParameters;
	//}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
