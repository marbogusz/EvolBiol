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
#include "heuristics/GuideTree.hpp"

using namespace std;

namespace EBC
{

class TripletSamplingTree
{
private:
	//all nodes
	vector<Node*> nodes;

	map<unsigned int, Node*> leafNodes;

	unordered_set<Node*> usedNodes;

	Node* mostRecentAncestor(Node* n1, Node* n2);

	double distanceBetween(Node* n1, Node* n2);

	DistanceMatrix* distMat;

	double idealTreeSize;

	double averageLeafbranch;
	double leafBranchSd;

	void fromNewick(const string& nString);

public:
	TripletSamplingTree(GuideTree& gt);

	//sample tripplets on a tree
	vector<array<unsigned int, 3> > sampleFromTree();

	//sample only based on the distance matrix
	vector<array<unsigned int, 3> > sampleFromDM();

	// node operator for comparision
};

} /* namespace EBC */

#endif /* TRIPLETSAMPLINGTREE_HPP_ */
