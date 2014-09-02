/*
 * TripletSamplingTree.hpp
 *
 *  Created on: Sep 1, 2014
 *      Author: root
 */

#ifndef TRIPLETSAMPLINGTREE_HPP_
#define TRIPLETSAMPLINGTREE_HPP_

#include <vector>
#include <unordered_set>
#include "Node.hpp"

using namespace std;

namespace EBC
{

class TripletSamplingTree
{
private:
	vector<Node*> nodes;

	vector<Node*> leafNodes;

	unordered_set<Node*> usedNodes;

	Node* mostRecentAncestor(Node* n1, Node* n2);

	double distanceBetween(Node* n1, Node* n2);



public:
	TripletSamplingTree();

	void fromNewick(string& nString);

	vector<Node*> sample();

	// node operator for comparision
};

} /* namespace EBC */

#endif /* TRIPLETSAMPLINGTREE_HPP_ */
