/*
 * ModelEstimator.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef MODELESTIMATOR_HPP_
#define MODELESTIMATOR_HPP_


#include <heuristics/PairSamplingTree.hpp>
#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include "core/ProgramException.hpp"
#include "core/PMatrixTriple.hpp"

#include "models/SubstitutionModelBase.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/GTRModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"

#include "heuristics/TripletAligner.hpp"
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

class ModelEstimator : public IOptimizable
{
protected:


	Dictionary* dict;

	SubstitutionModelBase* substModel;
	IndelModel* indelModel;

	Sequences* inputSequences;
	Maths* maths;

	PairSamplingTree pst;

	vector<ForwardPairHMM*> fwdHMMs;

	unsigned int gammaRateCategories;

	DistanceMatrix* distMat;

	void estimateParameters();

public:
	ModelEstimator(Sequences* inputSeqs, Definitions::ModelType model, DistanceMatrix* dm,
			unsigned int rateCategories, double alpha, bool estimateAlpha);

	virtual ~ModelEstimator();



	vector<double> getSubstitutionParameters();
	vector<double> getIndelParameters();
	double getAlpha();

	void recalculateHMMs();
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
