/*
 * GotohHMatrix.h
 *
 *  Created on: Feb 13, 2014
 *      Author: root
 */

#ifndef GOTOHHMATRIX_H_
#define GOTOHHMATRIX_H_

#include "heuristics/EvaluationMatrix.hpp"

namespace EBC
{
class GotohSMatrix;

class GotohHMatrix: public EBC::EvaluationMatrix
{
protected:

	virtual void initializeData();
	GotohSMatrix* sMatrix;

public:

	virtual double evaluate(unsigned int i, unsigned int j);

	GotohHMatrix(unsigned int xSize, unsigned int ySize, GotohScoringMatrix* sMat);

	void setSMatrix(GotohSMatrix* sm);
};

} /* namespace EBC */
#endif /* GOTOHHMATRIX_H_ */
