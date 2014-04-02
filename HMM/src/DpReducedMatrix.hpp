/*
 * DpReducedMatrix.h
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#ifndef DPREDUCEDMATRIX_H_
#define DPREDUCEDMATRIX_H_

#include <limits>
#include <iostream>

using namespace std;

namespace EBC
{

class DpReducedMatrix
{

private:

	void allocateData();

protected:

	unsigned int xSize, ySize;

	double minVal;
	double maxVal;

	virtual void initializeData() = 0;

public:

	//Two lines of data!
	double* buffer[2];

	//First and second row pointers
	double* previousRow;
	double* currentRow;

	DpReducedMatrix(unsigned int xSize, unsigned int ySize);

	virtual ~DpReducedMatrix();

	void setValue(unsigned int col, double value);

	double valueAtColumn(unsigned int col);

	double valueAtLeft(unsigned int col);

	double valueAtTop(unsigned int col);

	double valueAtDiagonal(unsigned int col);

	virtual void nextRow()=0;

};

} /* namespace EBC */
#endif /* DPREDUCEDMATRIX_H_ */
