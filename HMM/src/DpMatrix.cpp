/*
 * DpMatrix.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#include "DpMatrix.hpp"
#include <iostream>

using namespace std;

void EBC::DpMatrix::allocateData()
{
	matrixData = new TraceStep*[xSize];
	for(int i=0; i<xSize; i++)
	{
		matrixData[i] = new TraceStep[ySize];
	}
}

EBC::DpMatrix::DpMatrix(unsigned int xS, unsigned int yS) :
		xSize(xS), ySize(yS)
{
	maxVal = std::numeric_limits<double>::max();
	minVal = std::numeric_limits<double>::min();
	allocateData();
}

void EBC::DpMatrix::setValue(unsigned int x, unsigned int y, double value)
{
	//TODO - check bounds
	matrixData[x][y].score = value;
}


double EBC::DpMatrix::valueAt(unsigned int i, unsigned int j)
{
	return matrixData[i][j].score;
}

void EBC::DpMatrix::outputTrace()
{
	for(int i=0; i < xSize; i++)
	{
		for(int j=0; j < ySize; j++)
		{
			TraceStep& ts = matrixData[i][j];
					if (ts.vert &&  !ts.hor && !ts.diag)
						cout << "|";
					else if (ts.vert && ts.hor && !ts.diag)
						cout << "J";
					else if (ts.vert && ts.hor && ts.diag)
						cout << "*";
					else if (!ts.vert && ts.hor && !ts.diag)
						cout << "-";
					else if (!ts.vert && !ts.hor && ts.diag)
						cout << "\\";
					else if (ts.vert && ! ts.hor && ts.diag)
						cout << "V";
					else if (!ts.vert && ts.hor && ts.diag)
						cout << "7";
					else cout << ".";
		}
		cout << endl;
	}
}

void EBC::DpMatrix::setDiagonalAt(unsigned int i, unsigned int j)
{
	matrixData[i][j].diag = true;
}

void EBC::DpMatrix::setHorizontalAt(unsigned int i, unsigned int j)
{
	matrixData[i][j].hor = true;
}

void EBC::DpMatrix::setVerticalAt(unsigned int i, unsigned int j)
{
	matrixData[i][j].vert = true;
}
