/*
 * GuideTree.h
 *
 *  Created on: Aug 28, 2014
 *      Author: root
 */

#ifndef GUIDETREE_H_
#define GUIDETREE_H_


#include "core/Definitions.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include "core/HmmException.hpp"
#include "core/BioNJ.hpp"
#include <unordered_map>
#include <vector>
#include <string>

using namespace std;

namespace EBC
{

class GuideTree
{
protected:
	Dictionary* dict;
	Sequences* inputSequences;
	unsigned int kmerSize;
	unsigned int sequenceCount;
	vector<unordered_map<string,short>*>* kmers;
	vector<double> distances;

	vector<array<unsigned int, 3> > sampledTriplets;

	DistanceMatrix* distMat;

	string newickTree;


public:
	GuideTree(Sequences*);

	~GuideTree();

	void constructTree();

	DistanceMatrix* getDistanceMatrix()
	{
		return distMat;
	}

	const string& getNewickTree()
	{
		return newickTree;
	}

private:

	void extractKmers(string& seq, unordered_map<string,short>* umap);

	unsigned int commonKmerCount(unsigned int i, unsigned int j);

	double kimuraDist(double);

	inline double aaFunction(double kp)
	{
		//fitted to a function. Limit to 2 otherwise one might get weird results for distance matrix (negaive branch after NJ)
		double res = ((-0.0559362*kp)/(kp-0.9992776)) + ((-0.23315*kp)/(kp-1.3291755));
		if (std::isnan(res))
					res = 1.1;
		return min(res,1.1);
	}

	inline double nucFunction(double kP)
	{
		//double res = kP/pow(0.94015-3.14164*kP,0.3600898);
		//if (std::isnan(res))
		//	res = 1.1;
		//return min(res,1.1);
		return kP;
	}
};

} /* namespace EBC */

#endif /* GUIDETREE_H_ */
