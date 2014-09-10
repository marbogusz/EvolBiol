/*
 * Node.hpp
 *
 *  Created on: Sep 1, 2014
 *      Author: Marcin
 *
 *      Tree node class
 */

#include <vector>
#include <string>

#ifndef NODE_HPP_
#define NODE_HPP_

using namespace std;

namespace EBC
{

struct Node
{
public:

	unsigned int nodeId;

	Node* parent;

	bool leafNode;

	bool rootNode;

	unsigned int sequenceNo;

	string nodeName;

	double distanceToParent;

	friend bool operator== (Node &cP1, Node &cP2)
	{
		return cP1.nodeId == cP2.nodeId;
	}

	vector<Node*> children;

	Node(unsigned int id);

	void setName(string name);

	inline void setSequenceNumber(unsigned int no)
	{
		this->sequenceNo = no;
	}

	void setDistance(double distance);

	void setLeaf();

	void setParent(Node* pn);

	void setChild(Node* cn);

};

} /* namespace EBC */

#endif /* NODE_HPP_ */
