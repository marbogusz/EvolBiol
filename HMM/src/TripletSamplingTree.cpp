/*
 * TripletSamplingTree.cpp
 *
 *  Created on: Sep 1, 2014
 *      Author: root
 */

#include "TripletSamplingTree.hpp"
#include <regex>
#include <stack>
#include <random>

namespace EBC
{

Node* TripletSamplingTree::mostRecentAncestor(Node* n1, Node* n2)
{
	Node* tmpNode1;
	Node* tmpNode2;

	while(tmpNode1->parent != NULL)
	{
		tmpNode2 = n2;
		while(tmpNode2->parent != NULL)
		{
			if (tmpNode1 == tmpNode2)
				return tmpNode1;
			tmpNode2 = tmpNode2->parent;
		}
	tmpNode1 = tmpNode1->parent;
	}
	return NULL;
}

double TripletSamplingTree::distanceBetween(Node* n1, Node* n2)
{
	double distance = 0;
	Node* tmp = n1;

	if (n1 == n2)
		return 0;

	while (tmp != n2)
	{
		distance += tmp->distanceToParent;
		tmp = tmp->parent;
	}
	return distance;
}

TripletSamplingTree::TripletSamplingTree()
{
	// TODO Auto-generated constructor stub

}

void TripletSamplingTree::fromNewick(string& newick)
{
	//parse Newick!
	//start with a bracket, end with a bracket and a semicolon
	//TODO - check if the newick string is well formed
	bool endReached = false;
	unsigned int i = 0;
	unsigned int ids = 0;
	regex regDescDist("((\\w+):([0-9\\.]+)).*");
	regex regInterDist(":([0-9\\.]+).*");
	Node *tmpNode, *tmpParent, *tmpCurrent;
	stack<Node*> workNodes;

	string nodeName;
	double currentDistance;

	while (!endReached && i < newick.size())
	{
		if (std::isspace(newick[i]))
		{
			i++;
		}
		else if(newick[i] == '(')
		{
			tmpNode = new Node(++ids);
			nodes.push_back(tmpNode);
			workNodes.push(tmpNode);
			i++;
		}
		else if(newick[i] == ',')
		{
			workNodes.pop();
			if (workNodes.empty())
			{
				//push new father
				tmpNode = new Node(++ids);
				nodes.push_back(tmpNode);
				workNodes.push(tmpNode);
			}

			tmpParent = workNodes.top();
			tmpParent->setChild(tmpCurrent);
			tmpCurrent->setParent(tmpParent);

			tmpNode = new Node(++ids);
			nodes.push_back(tmpNode);
			workNodes.push(tmpNode);

			i++;
		}
		else if(newick[i] == ')')
		{

			tmpCurrent = workNodes.top();
			workNodes.pop();
			tmpParent = workNodes.top();
			tmpParent->setChild(tmpCurrent);
			tmpCurrent->setParent(tmpParent);
			//TODO - copy children!!

			i++;
		}
		else if(newick[i] == ';')
		{
			endReached = true;
			i++;
		}
		else
		{
			smatch basematch;
			string sstr = newick.substr(i);

			if (regex_match(sstr, basematch, regDescDist))
			{
				nodeName = basematch[2].str();
				currentDistance =  stof(basematch[3].str());
				tmpCurrent = workNodes.top();
				tmpCurrent->setName(nodeName);
				tmpCurrent->setDistance(currentDistance);
				tmpCurrent->setLeaf();
				this->leafNodes.pushBack(tmpCurrent);
				i += basematch[1].str().size();
			}
			else if (regex_match(sstr, basematch, regInterDist))
			{
				currentDistance =   stof(basematch[1].str());
				tmpCurrent = workNodes.top();
				tmpCurrent->setDistance(currentDistance);
				i += basematch[1].str().size() + 1;
			}
		}
	}
}

vector<Node*> TripletSamplingTree::sample()
{
	//randomly select 1 leaf
	double treeSize = 0;
	default_random_engine generator;
	uniform_int_distribution<int> distribution(0,leafNodes.size());
	distribution(generator);
	Node* first = leafNodes[distribution(generator)];
	//select the 2nd one from the distance matrix
	Node* second  = first;
	while (second == first)
	{
		second = leafNodes[distribution(generator)];
	}
	Node* ca12 =  mostRecentAncestor(first, second);

	treeSize = treeSize + distanceBetween(first, ca12) + distanceBetween(second, ca12);

	if (ca12->parent != NULL)
		this->usedNodes.emplace(ca12);

	//distance to root




	//get another one with distance between 0.3 and 0.7
	//what about root node - root node has no parent set - go through list and check ?????
	//find their common ancestor
	//based on the ancestor, calculate the distance to it and target the next taxon based on the distance to root + common to root
	//invalidate the leaves laying on the common path.

}

} /* namespace EBC */
