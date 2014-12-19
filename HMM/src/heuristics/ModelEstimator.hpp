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
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/GTRModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"

#include "heuristics/ModelEstimator.hpp"
#include "heuristics/GuideTree.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletAligner.hpp"
#include "heuristics/StateTransitionEstimator.hpp"
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

class ModelEstimator
{
protected:


	Dictionary* dict;

	SubstitutionModelBase* substModel;
	IndelModel* indelModel;

	Sequences* inputSequences;
	Maths* maths;

	GuideTree* gtree;
	TripletSamplingTree tst;

	StateTransitionEstimator* ste;
	SubstitutionModelEstimator* sme;

	TripletAligner* tal;
	ViterbiPairHMM* vphmm;

	vector<array<vector<SequenceElement>, 3> > tripleAlignments;
	vector<array<vector<SequenceElement>, 4> > pairAlignments;

	map<double, pair<vector<SequenceElement>, vector<SequenceElement> > > alSamples;

	double totalSampleLnl;

	vector<array<string, 3> > alignments;

	vector<array<unsigned int, 3> > tripletIdxs;

	unsigned int gammaRateCategories;

	bool estimateAlpha;

	void estimateTripleAlignment(Definitions::ModelType model, std::vector<double>& indel_params, std::vector<double>& subst_params, double al, double divergence);

	void sampleAlignments(ForwardPairHMM* hmm);

public:
	ModelEstimator(Sequences* inputSeqs, Definitions::ModelType model,
			Definitions::OptimizationType ot,
			unsigned int rateCategories, double alpha, bool estimateAlpha,std::vector<double> indel_params, std::vector<double> subst_params, double divergence);

	virtual ~ModelEstimator();

	vector<double> getSubstitutionParameters();
	vector<double> getIndelParameters();
	double getAlpha();

	GuideTree* getGuideTree()
	{
		return gtree;
	}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
