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

	unsigned int xBound;
	unsigned int yBound;

	unsigned int currentRowNo;

	inline void switchPointers()
	{
		double* tmp = firstRow;
		firstRow = secondRow;
		secondRow = tmp;
		currentRowNo++;
	}

protected:

	unsigned int xSize, ySize;

	double minVal;
	double maxVal;

	virtual void initializeData() = 0;

	virtual void setWholeRow(unsigned int row, double value);

	virtual void setWholeCol(unsigned int col, double value);


public:

	//Two lines of data!
	double* buffer[2];

	//First and second row pointers
	double* firstRow;
	double* secondRow;

	DpReducedMatrix(unsigned int xSize, unsigned int ySize);

	virtual ~DpReducedMatrix();

	void setValue(unsigned int x,unsigned int y, double value);

	double valueAt(unsigned int i, unsigned int j);

};

} /* namespace EBC */
#endif /* DPREDUCEDMATRIX_H_ */
