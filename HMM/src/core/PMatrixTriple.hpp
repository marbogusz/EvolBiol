/*
 * PMatrixTriple.hpp
 *
 *  Created on: Sep 17, 2014
 *      Author: root
 */

#ifndef PMATRIXTR_HPP_
#define PMATRIXTR_HPP_

#include "models/SubstitutionModelBase.hpp"
#include "core/HmmException.hpp"
#include "core/PMatrix.hpp"
#include <vector>
#include <array>

using namespace std;

namespace EBC
{

//TODO - split into 2 regular pt and pairwise
class PMatrixTriple : public PMatrix
{

public:
	PMatrixTriple(SubstitutionModelBase* m);
	virtual ~PMatrixTriple();


	void calculate();

	double getTransitionProb(unsigned int xi, unsigned int yi, unsigned int rateCat = 0);

	double getTripleSitePattern(unsigned int root,const array<unsigned char, 3>& nodes, PMatrixTriple* pm2, PMatrixTriple* pm3);

	void summarize();

};

} /* namespace EBC */

#endif /* PMATRIX_HPP_ */
