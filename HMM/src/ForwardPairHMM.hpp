/*
 * ForwardPairHMM.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef FORWARDPAIRHMM_HPP_
#define FORWARDPAIRHMM_HPP_

#include "EvolutionaryPairHMM.hpp"
#include <lbfgs.h>

namespace EBC
{

class ForwardPairHMM: public EBC::EvolutionaryPairHMM
{

private:

	//BFGS optimization wrapper for lbgfs
	class BFGS
	{
	protected:
		lbfgsfloatval_t* currentParameters;
		lbfgsfloatval_t likelihood;
		int paramsCount;
		double smallDiff;
		double smallDiffExp;
		double smallDiffExpInverse;
		unsigned int smallDiffInverse;

		double tempParams[100];// Params limit of 100;

		//lbfgsfloatval_t* parameterSpace;
		//parent pointer FIXME - pass only the params
		ForwardPairHMM* parent;

		static lbfgsfloatval_t _evaluate(void *instance, const lbfgsfloatval_t *x,
				lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
		{
			return reinterpret_cast<BFGS*>(instance)->calculateGradients(x, g, n, step);
		}

		static int _progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g,
				const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
				const lbfgsfloatval_t step,int n,int k, int ls)
		{
			return reinterpret_cast<BFGS*>(instance)->reportProgress(x, g, fx, xnorm, gnorm, step, n, k, ls);
		}

		lbfgsfloatval_t calculateGradients(const lbfgsfloatval_t *x,lbfgsfloatval_t *g,
				const int n, const lbfgsfloatval_t step);

		int reportProgress(const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
				const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls);

		void copyWithSmallDiff(int pos, const lbfgsfloatval_t* x);

		double logistic(double);

	public:
		BFGS(ForwardPairHMM* enclosing);
		virtual ~BFGS();
		void optimize();
	};


protected:

	virtual void initializeModels();

	BFGS* bfgs;

public:
	ForwardPairHMM(Sequences* inputSeqs, bool optimize);

	virtual ~ForwardPairHMM();

	double runForwardAlgorithm();

	double runForwardIteration(const double * bfgsParameters);
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
