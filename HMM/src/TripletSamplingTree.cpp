/*
 * TripletSamplingTree.cpp
 *
 *  Created on: Sep 1, 2014
 *      Author: root
 */

#include "TripletSamplingTree.hpp"
#include <regex>
#include <stack>

namespace EBC
{

Node* TripletSamplingTree::mostRecentAncestor(Node* n1, Node* n2)
{
}

double TripletSamplingTree::distanceBetween(Node* n1, Node* n2)
{
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
	unsigned int j = 0;
	regex regDescDist(R"((\w+):([0-9\.]))\)");
	regex regInterDist(R"([0-9\.])");
	Node* tmpNode, tmpParent;
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
			tmpNode = new Node();
			nodes.push_back(tmpNode);
			workNodes.push(tmpNode);
			i++;
		}
		else if(newick[i] == ',')
		{
			tmpNode = workNodes.top();
			workNodes.pop();

			i++;
		}
		else if(newick[i] == ')')
		{
			i++;
		}
		else
		{
			auto m = cmatch {};
			string sstr = newick.substr(i);
			if (regex_match(sstr, m, regDescDist))
			{
				//FIXME - length of the expression
				nodeName = "";
				currentDistance = 0;
				tmpNode = workNodes.top();
				workNodes.pop();
				tmpNode->setName(nodeName);
				tmpNode->setDistance(currentDistance);
				tmpParent = workNodes.top();
				tmpNode->setParent(tmpParent);
				tmpParent->addChild(tmpNode);

				i += 10;
			}
			else if (regex_match(sstr, m, regInterDist))
			{

			}
		}
	}



}

} /* namespace EBC */
