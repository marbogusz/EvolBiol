/*
 * MlEstimator.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef MLESTIMATOR_HPP_
#define MLESTIMATOR_HPP_


#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/MlEstimator.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "hmm/ViterbiPairHMM.hpp"
#include "models/SubstitutionModelBase.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "models/IndelModel.hpp"
#include "core/Sequences.hpp"
#include <sstream>
#include "core/HmmException.hpp"
#include "core/PMatrix.hpp"

#include <vector>
#include <map>
#include <array>

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
	//vector<SubstitutionModelBase*> substs;

	vector<PMatrix*> ptMatrices;

	OptimizedModelParameters* modelParams;

	vector<map<array<unsigned int, 2>, unsigned int> > patterns;

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
