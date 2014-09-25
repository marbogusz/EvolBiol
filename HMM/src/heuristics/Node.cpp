/*
 * Node.cpp
 *
 *  Created on: Sep 1, 2014
 *      Author: root
 */

#include "heuristics/Node.hpp"

namespace EBC
{

Node::Node(unsigned int id) : nodeId(id), parent(NULL), leafNode(false), rootNode(false), distanceToParent(0)
{
}

void Node::setName(string name)
{
	this->nodeName = name;
}

void Node::setDistance(double distance)
{
	this->distanceToParent = distance;
}

void Node::setLeaf()
{
	this->leafNode = true;
}

void Node::setParent(Node* pn)
{
	this->parent = pn;
}

void Node::setChild(Node* cn)
{
	this->children.insert(children.end(), cn->children.begin(), cn->children.end());
	this->children.push_back(cn);
}

} /* namespace EBC */
