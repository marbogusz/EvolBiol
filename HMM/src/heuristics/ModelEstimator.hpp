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

#include "heuristics/ModelEstimator.hpp"
#include "heuristics/GuideTree.hpp"
#include "heuristics/TripletSamplingTree.hpp"

#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletAligner.hpp"

#include "heuristics/StateTransitionEstimator.hpp"
#include "heuristics/SubstitutionModelEstimator.hpp"


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

	Sequences* inputSequences;
	Maths* maths;

	GuideTree gtree;
	TripletSamplingTree tst;

	StateTransitionEstimator* ste;
	SubstitutionModelEstimator* sme;

	vector<array<vector<SequenceElement>, 3> > tripleAlignments;

	vector<array<string, 3> > alignments;

	unsigned int gammaRateCategories;

	bool estimateAlpha;

public:
	ModelEstimator(Sequences* inputSeqs, Definitions::ModelType model,
			Definitions::OptimizationType ot,
			unsigned int rateCategories, double alpha, bool estimateAlpha);

	virtual ~ModelEstimator();

	vector<double> getSubstitutionParameters();
	vector<double> getIndelParameters();
	double getAlpha();

};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
