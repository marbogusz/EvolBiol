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
	buffer[0] = new double[xSize];
	buffer[1] = new double[xSize];
	firstRow = buffer;
	secondRow = buffer + 1;
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
	currentRowNo = 1;
}

void EBC::DpReducedMatrix::setValue(unsigned int x, unsigned int y, double value)
{
	//TODO - auto switch ?
	unsigned int diff = currentRowNo - x;
	switch(diff)
	{
	case 0:
		secondRow[y] = value;
		break;
	case 1:
		firstRow[y] = value;
	}
}

void EBC::DpReducedMatrix::setWholeRow(unsigned int row, double value)
{
	for (int i=0; i<ySize; i++)
	{

	}
}

void EBC::DpReducedMatrix::setWholeCol(unsigned int col, double value)
{
	for (int i=0; i<xSize; i++)
	{

	}
}


double EBC::DpReducedMatrix::valueAt(unsigned int i, unsigned int j)
{
	return secondRow[j];
}















