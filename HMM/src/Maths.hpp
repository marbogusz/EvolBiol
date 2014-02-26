/*
 * Maths.h
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef ALGEBRA_H_
#define ALGEBRA_H_

#include <ctime>
#include <cstdlib>
#include <cmath>


namespace EBC
{

//Various FUNCTIONS
class Maths
{
private:

	unsigned int z_rndu;

public:
	Maths();

	//Get a random within bounds
	double getRandom(double, double);

	//Exponentiate a reversible model matrix Q, store the result in P
	void revMatrixExponentiation(double* Q, double* P);

	//Shamelessly copied from Ziheng Yang's PAML
	int eigenRealSym(double A[], int n, double Root[], double work[]);

		//Shamelessly copied from Ziheng Yang's PAML
	int eigenQREV (double Q[], double pi[], int n,
			double Root[], double U[], double V[], double spacesqrtpi[]);

	void HouseholderRealSym(double a[], int n, double d[], double e[]);

	//Ziheng Yang's PAML fast random
	double rndu (void);

	void EigenSort(double d[], double U[], int n);

	int EigenTridagQLImplicit(double d[], double e[], int n, double z[]);

	//returns new matrix
	double* matrixMultiply(double *matA, double *matB, int size);

	//modifies an existing one
	void matrixByDiagonalMultiply(double *matA, double *matDiag, int size);

	//modifies existing one
	void vectorMultiply(double* vector, double factor, int size);

	void expLambdaT(double* lambda, double t, int size);

};

} /* namespace EBC */
#endif /* ALGEBRA_H_ */
