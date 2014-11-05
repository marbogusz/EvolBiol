/*
 * SubstitutionModelBase.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "models/SubstitutionModelBase.hpp"

namespace EBC
{

SubstitutionModelBase::SubstitutionModelBase(Dictionary* dict, Maths* alg, unsigned int rateCategories, unsigned int parameter_count)
	: dictionary(dict), maths(alg), rateCategories(rateCategories), paramsNumber(parameter_count),
	  parameterHiBounds(parameter_count), parameterLoBounds(parameter_count)
{
	this->matrixSize = dict->getAlphabetSize();
	this->matrixFullSize = matrixSize*matrixSize;
	this->allocateMatrices();
	this->meanRate=0.0;
	this-> roots = new double[matrixSize];
	//no alpha provided
	this->alpha = 0;
	//FIXME - maybe initial alpha set to a value ?
}


void SubstitutionModelBase::allocateMatrices()
{
	this->qMatrix = new double[matrixFullSize];
	this->uMatrix = new double[matrixFullSize];
	this->vMatrix = new double[matrixFullSize];
	this->squareRoots = new double[matrixFullSize];

	if (rateCategories > 1)
	{
		gammaFrequencies = new double[rateCategories];
		gammaRates = new double[rateCategories];
	}
	else
	{
		gammaFrequencies = new double[1];
		gammaRates = new double[1];
		gammaFrequencies[0] = gammaRates[0] = 1;
		meanRate = 1.0;
	}

	/*
	this->sitePatterns = new double*[matrixSize+1];
	this->siteProbabilities = new double*[matrixSize+1];
	for (int i =0; i<= matrixSize; i++ )
	{
		sitePatterns[i] = new double[matrixSize+1];
		siteProbabilities[i] = new double[matrixSize+1];;
	}
	*/
}


void SubstitutionModelBase::destroyMatrices()
{
	delete [] qMatrix;
	delete [] uMatrix;
	delete [] vMatrix;
	delete [] squareRoots;
	if (this->gammaFrequencies !=NULL )
		delete[] gammaFrequencies;
	if (this->gammaRates !=NULL )
		delete[] gammaRates;

}

void SubstitutionModelBase::doEigenDecomposition()
{
	std::fill(roots, roots+matrixSize, 0);
	std::fill(uMatrix, uMatrix+matrixFullSize, 0);
	std::fill(vMatrix, vMatrix+matrixFullSize, 0);
	std::fill(squareRoots, squareRoots+matrixFullSize, 0);
	this->maths->eigenQREV(qMatrix, piFreqs, matrixSize, roots, uMatrix, vMatrix, squareRoots);
}

double* SubstitutionModelBase::calculatePt(double t, unsigned int rateCategory)
{
	double *tmpRoots, *tmpUroots, *matrix;
	tmpRoots = maths->expLambdaT(roots, t*gammaRates[rateCategory], matrixSize);
	tmpUroots = maths->matrixByDiagonalMultiply(uMatrix, tmpRoots, matrixSize);
	matrix = this->maths->matrixMultiply(tmpUroots, vMatrix, matrixSize);
	delete[] tmpRoots;
	delete[] tmpUroots;
	return matrix;
}

void SubstitutionModelBase::setDiagonalMeans()
{
		int i,j;
		double sum;

		maths->matrixByDiagonalMultiplyMutable(qMatrix,piFreqs,matrixSize);
		meanRate = 0;

		for (i=0; i< this->matrixSize; i++)
		{
			qMatrix[i*matrixSize+i] = 0;
			sum = 0.0;
			for (j=0; j < this->matrixSize; j++)
			{
				sum -= qMatrix[i*matrixSize + j];
			}
			qMatrix[i*matrixSize + i] = sum;
			meanRate -= sum*piFreqs[i];
		}
		for (i=0; i< this->matrixSize; i++)
		{
				for (j=0; j < this->matrixSize; j++)
				{
					qMatrix[i*matrixSize + j] /= meanRate;
				}
		}

}

SubstitutionModelBase::~SubstitutionModelBase()
{
	delete[] this-> roots;
	destroyMatrices();

	//site patterns
	/*
	for (int i =0; i<= matrixSize; i++ )
	{
			delete[] sitePatterns[i];
			delete[] siteProbabilities[i];
	}
	delete[] sitePatterns;
	delete[] siteProbabilities;
	 */
}

void SubstitutionModelBase::calculateGamma()
{
	int useMedian = 0;
	if(!(alpha <=0 || rateCategories == 1))
	{
		this->maths->DiscreteGamma(gammaFrequencies, gammaRates, alpha, alpha, rateCategories, useMedian);
	}
}

void SubstitutionModelBase::setObservedFrequencies(double* observedFrequencies)
{
	//Size and order is clear
	this->piFreqs = observedFrequencies;

}

double SubstitutionModelBase::getEquilibriumFrequencies(unsigned int xi)
{
	return piFreqs[xi];
}

void SubstitutionModelBase::summarizeRates()
{
	if(this->rateCategories != 1)
	{
		cout << "Rate categories: " << endl;
		for(int r=0; r< rateCategories; r++)
		{
			cout << gammaRates[r] << "\t\t";
		}
		cout << endl;
		cout << "Rate frequencies : " << endl;
		for(int f=0; f< rateCategories; f++)
		{
			cout << this->gammaFrequencies[f] << "\t\t";
		}
		cout << endl << "Alpha : " << alpha;
		cout << endl;
	}

	/*
	cout << "Rate Matrix: " << endl;

		for (int i =0; i< matrixSize; i++ )
		{
				for (int j =0; j<  matrixSize; j++ )
				{
					cout << qMatrix[(i*matrixSize)+j] << "\t\t";
				}
				cout << endl;
		}
		cout << endl;
*/
}

/*
 *
 void SubstitutionModelBase::calculateGammaPtMatrices()
{
	double *tmpRoots, *tmpUroots;
	unsigned int i;
	int matrixCount = rateCategories == 1 ? 1 : rateCategories;

	std::vector<double*> pMatVector(matrixCount);

	for(i = 0; i< matrixCount; i++)
	{

		tmpRoots = maths->expLambdaT(roots, time*gammaRates[i], matrixSize);
		tmpUroots = maths->matrixByDiagonalMultiply(uMatrix, tmpRoots, matrixSize);
		//tmpPmatrix = this->maths->matrixMultiply(tmpUroots, vMatrix, matrixSize);
		pMatVector[i] = this->maths->matrixMultiply(tmpUroots, vMatrix, matrixSize);
		delete[] tmpRoots;
		delete[] tmpUroots;
	}
	//now we have a vector of p-mats, combine into 1 matrix since the probs are equal in this model
	for (i=1; i<matrixCount; i++)
	{
		maths->matrixAppend(pMatVector[0],pMatVector[i],matrixSize);
		delete[] pMatVector[i];
	}
	this->pMatrix = pMatVector[0];
	if(rateCategories != 1)
	{
		this->maths->matrixScale(pMatrix,gammaFrequencies[0],matrixSize);
	}
}

void SubstitutionModelBase::calculateSitePatterns()
{
	//includes gaps - does not discard missing data!
	for (int i =0; i<= matrixSize; i++ )
		for (int j =0; j<= matrixSize; j++ )
		{
			if (i == matrixSize)
			{
				siteProbabilities[i][j] = this->piFreqs[j];
				sitePatterns[i][j]  = log(this->piFreqs[j]);
			}
			else if (j == matrixSize)
			{
				siteProbabilities[i][j] = this->piFreqs[i];
				sitePatterns[i][j]  = log(this->piFreqs[i]);
			}
			else
			{
				siteProbabilities[i][j] = this->getPiXiPXiYi(i,j);
				sitePatterns[i][j] = log(this->getPiXiPXiYi(i,j));

			}
	}
	sitePatterns[matrixSize][matrixSize] = 0;
	siteProbabilities[matrixSize][matrixSize] = 0;
}

void SubstitutionModelBase::summarizePatterns()
{
	if(this->rateCategories != 0)
	{
		cout << "Rate categories: " << endl;
		for(int r=0; r< rateCategories; r++)
		{
			cout << gammaRates[r] << "\t\t";
		}
		cout << endl;
	}

	cout << "Rate Matrix: " << endl;

		for (int i =0; i< matrixSize; i++ )
		{
				for (int j =0; j<  matrixSize; j++ )
				{
					cout << qMatrix[(i*matrixSize)+j] << "\t\t";
				}
				cout << endl;
		}
		cout << endl;


	cout << "Site patterns " << endl;
	for (int i =0; i<= matrixSize; i++ )
	{
			for (int j =0; j<= matrixSize; j++ )
			{
				cout << sitePatterns[i][j] << "\t\t";
			}
			cout << endl;
	}

	cout << endl <<"Site probabilities " << endl;
	for (int i =0; i<= matrixSize; i++ )
	{
			for (int j =0; j<= matrixSize; j++ )
			{
				cout << siteProbabilities[i][j] << "\t\t";
			}
			cout << endl;
	}
	cout << endl << "P(t) matrix" << endl;
	for (int i =0; i< matrixSize; i++ )
	{
			for (int j =0; j<  matrixSize; j++ )
			{
				cout << pMatrix[(i*matrixSize)+j] << "\t\t";
			}
			cout << endl;
	}
	cout << endl;


}

double SubstitutionModelBase::getPiXiPXiYi(unsigned int xi, unsigned int yi)
{
	//FIXME - RATES
		return piFreqs[xi]*pMatrix[(xi*matrixSize)+yi];
}

double SubstitutionModelBase::getPXiYi(unsigned int xi, unsigned int yi)
{
	//FIXME - RATES
		return pMatrix[(xi*matrixSize)+yi];
}



double SubstitutionModelBase::getSitePattern(unsigned int xi, unsigned int yi)
{
	return this->sitePatterns[xi][yi];
}

double SubstitutionModelBase::getSiteProbability(unsigned int xi, unsigned int yi)
{
	return this->siteProbabilities[xi][yi];
}
*/

} /* namespace EBC */


