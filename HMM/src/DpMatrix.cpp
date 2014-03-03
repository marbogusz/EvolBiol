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

EBC::DpMatrix::~DpMatrix()
{
	for(int i=0; i<xSize; i++)
	{
		delete[] matrixData[i];
	}
	delete matrixData;
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

void EBC::DpMatrix::setWholeRow(unsigned int row, double value)
{
	for (int i=0; i<ySize; i++)
	{
		matrixData[row][i].score = value;
	}
}

void EBC::DpMatrix::setWholeCol(unsigned int col, double value)
{
	for (int i=0; i<xSize; i++)
	{
		matrixData[i][col].score = value;
	}
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

void EBC::DpMatrix::traceback(string& seq_a, string& seq_b, std::pair<string,string>* alignment)
{
	unsigned int i = xSize-1;
	unsigned int j = ySize-1;

	while(i>0 || j >0)
	{
		if (matrixData[i][j].diag)
		{
			alignment->first += seq_a[i-1];
			alignment->second += seq_b[j-1];
			i--;
			j--;
		}
		else if (matrixData[i][j].hor)
		{
			alignment->second += seq_b[j-1];
			alignment->first += '-';
			j--;
		}
		//vert
		else
		{
			alignment->first += seq_a[i-1];
			alignment->second += '-';
			i--;
		}
	}
}













