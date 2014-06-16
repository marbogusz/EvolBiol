/*
 * DpMatrixLoMem.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#include "DpMatrixLoMem.hpp"
#include <iostream>

using namespace std;

void EBC::DpMatrixLoMem::allocateData()
{

	buffer[0] = new double[ySize];
	buffer[1] = new double[ySize];

	previousRow = buffer[0];
	currentRow = buffer[1];

	for (unsigned int i=0; i< ySize; i++)
	{
		currentRow[i] = previousRow[i] = minVal;
	}
	currentRowIndex = 0;
}

EBC::DpMatrixLoMem::~DpMatrixLoMem()
{
	delete[] buffer[0];
	delete[] buffer[1];
}

EBC::DpMatrixLoMem::DpMatrixLoMem(unsigned int xS, unsigned int yS) :
		xSize(xS), ySize(yS)
{
	maxVal = std::numeric_limits<double>::max();
	minVal = -10000;//std::numeric_limits<double>::min();
	allocateData();
	//set to the beginning;
	//always associated with the second ptr.
}

void EBC::DpMatrixLoMem::setValue(unsigned int x,unsigned int y, double value)

double EBC::DpMatrixLoMem::valueAt(unsigned int i, unsigned int j)
{
	if(i==currentRowIndex)
	{
		return currentRow[j];
	}
	else if((i==currentRowIndex-1))
	{
		return previousRow[j];
	}
}

void EBC::DpMatrixLoMem::setSrc(unsigned int i, unsigned int j, DpMatrixBase*){}

void EBC::DpMatrixLoMem::setDiagonalAt(unsigned int i, unsigned int j){}

void EBC::DpMatrixLoMem::setHorizontalAt(unsigned int i, unsigned int j){}

void EBC::DpMatrixLoMem::setVerticalAt(unsigned int i, unsigned int j){}

void EBC::DpMatrixLoMem::setWholeRow(unsigned int row, double value)
{
	if(row == currentRowIndex)
	{
		std::fill(currentRow,currentRow+ySize, value);
		currentRowIndex++;
	}
}



void EBC::DpMatrixLoMem::setWholeCol(unsigned int col, double value)
{
	this->buffer[0] = this->buffer[1] = minVal;
}

void EBC::DpMatrixLoMem::nextRow()
{
	double* tmp = previousRow;
	previousRow = currentRow;
	currentRow = tmp;
	std::fill(currentRow,currentRow+ySize, -10000.0);
}

void EBC::DpMatrixLoMem::setValue(unsigned int col, double value)
{
		currentRow[col] = value;
}

double EBC::DpMatrixLoMem::valueAtColumn(unsigned int col)
{
	return currentRow[col];
}

double EBC::DpMatrixLoMem::valueAtLeft(unsigned int col)
{
	return currentRow[col-1];
}

double EBC::DpMatrixLoMem::valueAtTop(unsigned int col)
{
	return previousRow[col];
}

double EBC::DpMatrixLoMem::valueAtDiagonal(unsigned int col)
{
	return previousRow[col-1];
}

















