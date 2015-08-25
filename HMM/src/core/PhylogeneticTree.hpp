/*
 * PhylogeneticTree.hpp
 *
 *  Created on: Sep 1, 2014
 *      Author: root
 */

#ifndef TREE_HPP_
#define TREE_HPP_

#include <vector>
#include <array>
#include <map>
#include <unordered_set>
#include "heuristics/Node.hpp"

using namespace std;

namespace EBC
{

class PhylogeneticTree
{
private:
	//all nodes
	vector<Node*> nodes;

	map<unsigned int, Node*> leafNodes;

	map<unsigned int, Node*> availableNodes;

	unordered_set<Node*> usedNodes;

	Node* mostRecentAncestor(Node* n1, Node* n2);

	double distanceToParent(Node* n1, Node* par);

public:
	PhylogeneticTree();

	~PhylogeneticTree();

	void fromNewick(const string& nString);

	double distanceById(unsigned int n1, unsigned int n2);
	double distanceByName(string& n1, string& n2);

};

} /* namespace EBC */

#endif /* TRIPLETSAMPLINGTREE_HPP_ */
