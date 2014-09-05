/*
 * GotohSMatrix.h
 *
 *  Created on: Feb 13, 2014
 *      Author: root
 */

#ifndef GOTOHSMATRIX_H_
#define GOTOHSMATRIX_H_

#include "heuristics/EvaluationMatrix.hpp"


namespace EBC
{

//forwarded declarations
class GotohVMatrix;
class GotohHMatrix;

class GotohSMatrix: public EBC::EvaluationMatrix
{

protected:

	virtual void initializeData();
	GotohVMatrix* vMatrix;
	GotohHMatrix* hMatrix;

public:

	virtual double evaluate(unsigned int i, unsigned int j);

	GotohSMatrix(unsigned int xSize, unsigned int ySize, GotohScoringMatrix* sMat);

	void setHMatrix(GotohHMatrix* sm);

	void setVMatrix(GotohVMatrix* vm);

	void traceback();


};

} /* namespace EBC */
#endif /* GOTOHSMATRIX_H_ */
