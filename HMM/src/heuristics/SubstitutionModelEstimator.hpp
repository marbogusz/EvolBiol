/*
 * SubstitutionModelEstimator.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef SUBSTMODELESTIMATOR_HPP_
#define SUBSTMODELESTIMATOR_HPP_


#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/HmmException.hpp"
#include "core/PMatrixTriple.hpp"
#include "core/IOptimizable.hpp"
#include "core/Optimizer.hpp"
#include "core/Sequences.hpp"

#include "models/SubstitutionModelBase.hpp"

#include <vector>
#include <array>

using namespace std;


namespace EBC
{

class SubstitutionModelEstimator : public IOptimizable
{

protected:

	Optimizer* bfgs;

	Dictionary* dict;

	Sequences* inputSequences;

	Maths* maths;

	SubstitutionModelBase* substModel;

	vector<array<PMatrixTriple* ,3> > ptMatrices;

	vector<map<array<unsigned char, 3>, unsigned int> > patterns;

	unsigned int gammaRateCategories;

	OptimizedModelParameters* modelParams;

	double alpha;

	bool estimateSubstitutionParams;
	bool estimateAlpha;

	unsigned int currentTriplet;

public:
	SubstitutionModelEstimator(Sequences* inputSeqs, SubstitutionModelBase* model,
			Definitions::OptimizationType ot,unsigned int rateCategories, double alpha,
			bool estimateAlpha, unsigned int matCount);

	virtual ~SubstitutionModelEstimator();

	void addTriplet(array<vector<unsigned char>*, 3> tripleAlignment, unsigned int tiplet);

	double runIteration();

	void optimize();

	void clean();

	double getTripletDivergence(unsigned int triplet, unsigned int branch)
	{
		//no bound checks - beware
		return modelParams->getDivergenceTime(((triplet*3)+branch));
	}

	OptimizedModelParameters* getModelParams()
	{
		return modelParams;
	}

	SubstitutionModelBase* getSubstModel()
	{
		return substModel;
	}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
