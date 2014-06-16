/*
 * DpReducedMatrix.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#include "DpReducedMatrix.hpp"
#include <iostream>

using namespace std;

void EBC::DpReducedMatrix::allocateData()
{

	buffer[0] = new double[ySize];
	buffer[1] = new double[ySize];
	previousRow = buffer[0];
	currentRow = buffer[1];

}

EBC::DpReducedMatrix::~DpReducedMatrix()
{
	delete[] buffer[0];
	delete[] buffer[1];
}

EBC::DpReducedMatrix::DpReducedMatrix(unsigned int xS, unsigned int yS) :
		xSize(xS), ySize(yS)
{
	maxVal = std::numeric_limits<double>::max();
	minVal = -10000;//std::numeric_limits<double>::min();
	allocateData();
	//set to the beginning;
	//always associated with the second ptr.
}

void EBC::DpReducedMatrix::setValue(unsigned int col, double value)
{
		currentRow[col] = value;
}

double EBC::DpReducedMatrix::valueAtColumn(unsigned int col)
{
	return currentRow[col];
}

double EBC::DpReducedMatrix::valueAtLeft(unsigned int col)
{
	return currentRow[col-1];
}

double EBC::DpReducedMatrix::valueAtTop(unsigned int col)
{
	return previousRow[col];
}

double EBC::DpReducedMatrix::valueAtDiagonal(unsigned int col)
{
	return previousRow[col-1];
}















