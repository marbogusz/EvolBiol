/*
 * ModelEstimator.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef RAWINDELMODELESTIMATOR_HPP_
#define RAWINDELMODELESTIMATOR_HPP_


#include <heuristics/PairSamplingTree.hpp>
#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include "core/ProgramException.hpp"
#include "core/PMatrixTriple.hpp"
#include "core/IOptimizable.hpp"

#include "models/SubstitutionModelBase.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/GTRModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"

#include "heuristics/TripletAligner.hpp"
#include "heuristics/SubstitutionModelEstimator.hpp"
#include "heuristics/PairSamplingTree.hpp"

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

class RawIndelEstimator : public IOptimizable
{
protected:

	Optimizer* bfgs;
	OptimizedModelParameters* modelParams;

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

	vector<array<unsigned int, 2> > pairIdxs;
	vector<double> pairwiseTimes;
	vector<Band*> bands;


	unsigned int pairIdxsSize;

	double maxTime;



public:
	RawIndelEstimator(Sequences* inputSeqs, Definitions::ModelType model, DistanceMatrix* dm);

	virtual ~RawIndelEstimator();

	vector<double> getIndelParameters();
	double getAlpha();

	double runIteration();

	void calculateInitials(Definitions::ModelType model);

	OptimizedModelParameters* getModelParams()
	{
		return modelParams;
	}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
