/*
 * GotohVMatrix.h
 *
 *  Created on: Feb 13, 2014
 *      Author: root
 */

#ifndef GOTOHVMATRIX_H_
#define GOTOHVMATRIX_H_

#include "heuristics/EvaluationMatrix.hpp"

namespace EBC
{
class GotohSMatrix;

class GotohVMatrix: public EBC::EvaluationMatrix
{
protected:

	virtual void initializeData();
	GotohSMatrix* sMatrix;

public:

	virtual double evaluate(unsigned int i, unsigned int j);

	GotohVMatrix(unsigned int xSize, unsigned int ySize, GotohScoringMatrix* sMat);

	void setSMatrix(GotohSMatrix* SM);

};

} /* namespace EBC */
#endif /* GOTOHVMATRIX_H_ */
