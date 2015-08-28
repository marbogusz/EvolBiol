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
#include <vector>

#include "core/SequenceElement.hpp"
#include "core/Dictionary.hpp"
#include "core/Definitions.hpp"

using namespace std;

namespace EBC
{

class DpMatrixBase
{

protected:

	unsigned int xSize, ySize;

	double minVal;

	virtual void allocateData()=0;
public:

	virtual void setWholeRow(unsigned int row, double value)=0;

	virtual void setWholeCol(unsigned int col, double value)=0;

	DpMatrixBase(unsigned int xSize, unsigned int ySize)
	{
		this->xSize = xSize;
		this->ySize = ySize;
		this->minVal = Definitions::minMatrixLikelihood;
	}

	virtual ~DpMatrixBase() {}

	virtual void setValue(unsigned int x,unsigned int y, double value)=0;

	virtual double valueAt(unsigned int i, unsigned int j)=0;

	virtual void setSrc(unsigned int i, unsigned int j, DpMatrixBase*)=0;

	virtual void setDiagonalAt(unsigned int i, unsigned int j)=0;

	virtual void setHorizontalAt(unsigned int i, unsigned int j)=0;

	virtual void setVerticalAt(unsigned int i, unsigned int j)=0;

	virtual void traceback(string& seq_a, string& seq_b, std::pair<string,string>* alignment)=0;

	virtual void tracebackRaw(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dict, vector<std::pair<unsigned int, unsigned int> >&)=0;
};

} /* namespace EBC */
#endif /* DPMATRIXBASE_H_ */
