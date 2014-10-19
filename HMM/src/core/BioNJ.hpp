/*
 * BioNJ.hpp
 *
 *  Created on: Jun 24, 2014
 *      Author: root
 */

#ifndef BIONJ_H_
#define BIONJ_H_

#include <cstdio>
#include <cstdlib>
#include <string>
#include <ctime>
#include <vector>
#include "DistanceMatrix.hpp"


#define PREC 8                             /* precision of branch-lengths  */
#define PRC  100
#define LEN  1000                            /* length of taxon names        */

using namespace std;

typedef struct word
{
  char name[LEN];
  struct word *suiv;
}WORD;

typedef struct pointers
{
  WORD *head;
  WORD *tail;
}POINTERS;

namespace EBC
{

class BioNJ
{

private:
	unsigned int taxas;
	unsigned int pairs;
	vector<string> names;
	DistanceMatrix* times;
	vector<double> timesVec;
	double treeLength;

public:
	BioNJ(unsigned int size, DistanceMatrix* divergenceTimes);
	
	BioNJ(unsigned int size, vector<double> divergenceTimes);

	double getDist(unsigned int i, unsigned int j)
	{
		if (times != NULL)
			return times->getDistance(i,j);
		
		unsigned int a,b;
		unsigned int index;
		a=i;
		b=j;
		if (i>j)
		{
			b=i;
			a=j;
		}
		if (i==j) return 0;
		else
		{
			index = (b - a - 1) + (a*taxas) - (((1+a)/2.0)*(a*1.0));
			return timesVec[index];
		}
		
	}

	void   Initialize(float **delta, int n, POINTERS *trees);

	void   Compute_sums_Sx(float **delta, int n);

	void   Best_pair(float **delta, int r, int *a, int *b, int n);

	void   Finish(float **delta, int n, POINTERS *trees, stringstream& output);

	void   Concatenate(char chain1[LEN], int ind, POINTERS *trees, int post);

	void   Print_output(int i, POINTERS *trees, stringstream& output);

	float Distance(int i, int j, float **delta);

	float Variance(int i, int j, float **delta);

	float Sum_S(int i, float **delta);

	float Agglomerative_criterion(int i, int j, float **delta, int r);

	float Branch_length(int a, int b, float **delta, int r);

	float Reduction4(int a, float la, int b, float lb, int i, float lamda,
			 float **delta);

	float Reduction10(int a, int b, int i, float lamda, float vab, float
			  **delta);
	float Lamda(int a, int b, float vab, float **delta, int n, int r);

	float Finish_branch_length(int i, int j, int k, float **delta);

	int    Emptied(int i, float **delta);

	int    Symmetrize(float **delta, int n);

	string calculate();
};

} /* namespace EBC */

#endif /* BIONJ_H_ */












