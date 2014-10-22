/*
 * DpMatrixLoMem.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#include "hmm/DpMatrixLoMem.hpp"
#include <iostream>
#include "core/HmmException.hpp"

using namespace std;

void EBC::DpMatrixLoMem::allocateData()
{

	buffer[0] = new double[ySize];
	buffer[1] = new double[ySize];

	clear();

	for (unsigned int i=0; i< ySize; i++)
	{
		currentRow[i] = previousRow[i] = minVal;
	}
}

void EBC::DpMatrixLoMem::clear()
{
	previousRow = buffer[0];
	currentRow = buffer[1];
	currentRowIndex = 0;
	nextRowIndex =1;
	previousRowIndex = -1;
}

EBC::DpMatrixLoMem::~DpMatrixLoMem()
{
	delete[] buffer[0];
	delete[] buffer[1];
}

EBC::DpMatrixLoMem::DpMatrixLoMem(unsigned int xS, unsigned int yS) : DpMatrixBase(xS,yS)
{
	this->allocateData();
}

void EBC::DpMatrixLoMem::setValue(unsigned int i,unsigned int j, double value)
{
	if(i==nextRowIndex)
			nextRow();

	if(i==currentRowIndex)
	{
		currentRow[j] = value;
	}
	else if(i==previousRowIndex)
	{
		previousRow[j] = value;
	}

	else
	{
		string msg = "ERROR setValue() index out of bounds, x: " + std::to_string(i) + " y : " + std::to_string(j) + " current row is " + std::to_string(currentRowIndex) + "\n";
		throw HmmException(msg);
	}
}

double EBC::DpMatrixLoMem::valueAt(unsigned int i, unsigned int j)
{
	if(i==nextRowIndex)
		nextRow();

	if(i==currentRowIndex)
	{
		return currentRow[j];
	}
	//FIXME - check casting to int
	else if(i==previousRowIndex)
	{
		return previousRow[j];
	}
	else
	{
		string msg = "Value At : matrix index out of bounds, x: " + std::to_string(i) + " y : " + std::to_string(j) + " current row is " + std::to_string(currentRowIndex) + "\n";
		throw HmmException(msg);
	}
}

void EBC::DpMatrixLoMem::setSrc(unsigned int i, unsigned int j, DpMatrixBase*){}

void EBC::DpMatrixLoMem::setDiagonalAt(unsigned int i, unsigned int j){}

void EBC::DpMatrixLoMem::setHorizontalAt(unsigned int i, unsigned int j){}

void EBC::DpMatrixLoMem::setVerticalAt(unsigned int i, unsigned int j){}

void EBC::DpMatrixLoMem::setWholeRow(unsigned int row, double value)
{
	clear();
	if(row == currentRowIndex)
	{
		std::fill(currentRow,currentRow+ySize, value);
	}
}

void EBC::DpMatrixLoMem::setWholeCol(unsigned int col, double value)
{
	clear();
	this->buffer[0][0] = this->buffer[1][0] = this->minVal;
}

void EBC::DpMatrixLoMem::nextRow()
{
	double* tmp = previousRow;
	previousRow = currentRow;
	currentRow = tmp;
	std::fill(currentRow,currentRow+ySize, -10000.0);
	currentRowIndex++;
	nextRowIndex++;
	previousRowIndex++;
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

















