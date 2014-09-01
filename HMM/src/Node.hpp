/*
 * Node.hpp
 *
 *  Created on: Sep 1, 2014
 *      Author: Marcin
 *
 *      Tree node class
 */

#include <list>
#include <string>

#ifndef NODE_HPP_
#define NODE_HPP_

using namespace std;

namespace EBC
{

class Node
{
private:
	unsigned int nodeId;

	Node* parent;

	bool leafNode;

	bool rootNode;

	string* sequence;

	string nodeName;

	double distanceToParent;

	list<Node*> children;

public:
	Node();
};

} /* namespace EBC */

#endif /* NODE_HPP_ */
