/*
 * TripletSamplingTree.hpp
 *
 *  Created on: Sep 1, 2014
 *      Author: root
 */

#ifndef TRIPLETSAMPLINGTREE_HPP_
#define TRIPLETSAMPLINGTREE_HPP_

#include <vector>
#include <array>
#include <map>
#include <unordered_set>
#include "heuristics/Node.hpp"
#include "core/DistanceMatrix.hpp"

using namespace std;

namespace EBC
{

class PairSamplingTree
{
private:
	//all nodes
	vector<Node*> nodes;

	map<unsigned int, Node*> leafNodes;

	map<unsigned int, Node*> availableNodes;

	unordered_set<Node*> usedNodes;

	Node* mostRecentAncestor(Node* n1, Node* n2);

	double distanceToParent(Node* n1, Node* par);

	inline bool isWithinRange(double val, double ran)
	{
		if(val <= ran*1.2 || val >= ran*0.8)
			return true;
		return false;
	}

	inline bool isWithinRange(double val, pair<double, double> ran)
	{
		if(val >= ran.first && val <= ran.second)
			return true;
		return false;
	}

	DistanceMatrix* distMat;

	double idealTreeSize;

	double averageLeafbranch;
	double leafBranchSd;

public:
	PairSamplingTree(DistanceMatrix* dm);

	~PairSamplingTree();

	//sample only based on the distance matrix
	vector<array<unsigned int, 2> > sampleFromTree();

	// node operator for comparision
};

} /* namespace EBC */

#endif /* TRIPLETSAMPLINGTREE_HPP_ */
