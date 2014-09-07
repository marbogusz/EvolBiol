/*
 * TripletModelEstimator.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef TRIPLETMODELESTIMATOR_HPP_
#define TRIPLETMODELESTIMATOR_HPP_


#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "heuristics/TripletModelEstimator.hpp"
#include "models/SubstitutionModelBase.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include <sstream>
#include "core/HmmException.hpp"

#include <vector>

using namespace std;


namespace EBC
{

class TripletModelEstimator
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

		TripletModelEstimator* parent;

		Definitions::OptimizationType optimizationType;

	public:
		BFGS(TripletModelEstimator* enclosing, Definitions::OptimizationType ot);
		virtual ~BFGS();
		void optimize();

		double objectiveFunction(const column_vector& m);

		const column_vector objectiveFunctionDerivative(const column_vector& m);
	};


protected:

	BFGS* bfgs;

	Dictionary* dict;
	SubstitutionModelBase* substModel;
	Sequences* inputSequences;
	Maths* maths;

	unsigned int gammaRateCategories;

	bool estimateSubstitutionParams;
	bool estimateDivergence;
	bool estimateAlpha;

	OptimizedModelParameters* modelParams;

public:
	TripletModelEstimator(Sequences* inputSeqs, Definitions::ModelType model,
			Definitions::OptimizationType ot,
			unsigned int rateCategories, double alpha, bool estimateAlpha);

	virtual ~TripletModelEstimator();

	void runIteration(const column_vector& bfgsParameters);

	double runIteration();
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
