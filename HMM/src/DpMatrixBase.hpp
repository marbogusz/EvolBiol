/*
 * DpMatrixBase.h
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#ifndef DPMATRIXBASE_H_
#define DPMATRIXBASE_H_

#include <limits>
#include <iostream>

using namespace std;

namespace EBC
{

class DpMatrixBase
{

protected:

	unsigned int xSize, ySize;

	double minVal;
	double maxVal;

public:

	virtual void setWholeRow(unsigned int row, double value)=0;

	virtual void setWholeCol(unsigned int col, double value)=0;

	DpMatrixBase(unsigned int xSize, unsigned int ySize);

	virtual ~DpMatrixBase();

	virtual void setValue(unsigned int x,unsigned int y, double value)=0;

	virtual double valueAt(unsigned int i, unsigned int j)=0;

	virtual void setSrc(unsigned int i, unsigned int j, DpMatrixBase*)=0;

	virtual void setDiagonalAt(unsigned int i, unsigned int j)=0;

	virtual void setHorizontalAt(unsigned int i, unsigned int j)=0;

	virtual void setVerticalAt(unsigned int i, unsigned int j)=0;
};

} /* namespace EBC */
#endif /* DPMATRIXBASE_H_ */
