/*
 * PhylogeneticTree.cpp
 *
 *  Created on: Sep 1, 2014
 *      Author: root
 */

#include "core/PhylogeneticTree.hpp"
#include <regex>
#include <stack>
#include <random>
#include <iostream>
#include <cmath>
#include "core/Definitions.hpp"

namespace EBC
{

PhylogeneticTree::PhylogeneticTree(Sequences* seqs) : inputSeqs(seqs), dm(NULL)
{
	DEBUG("Creating PhylogeneticTree");

}

PhylogeneticTree::~PhylogeneticTree()
{
	for (auto nod : nodes)
		delete nod;
}


void PhylogeneticTree::buildDM()
{
	int sequenceCount = inputSeqs->getSequenceCount();
	dm = new DistanceMatrix(sequenceCount);
	for(int i = 0; i< sequenceCount; i++)
		for(int j = i+1; j< sequenceCount; j++)
		{
			dm->addDistance(i,j,distanceById(i,j));
		}
}

Node* PhylogeneticTree::mostRecentAncestor(Node* n1, Node* n2)
{
	Node* tmpNode1 = n1;
	Node* tmpNode2 = n2;

	while(tmpNode1 != NULL)
	{
		tmpNode2 = n2;
		while(tmpNode2 != NULL)
		{
			if (*tmpNode1 == *tmpNode2)
				return tmpNode1;
			tmpNode2 = tmpNode2->parent;
		}
	tmpNode1 = tmpNode1->parent;
	}
	return NULL;
}

double PhylogeneticTree::distanceToParent(Node* n1, Node* par)
{
	double distance = 0;
	Node* tmp = n1;

	if (*n1 == *par)
		return 0;

	while (tmp != par)
	{
		distance += tmp->distanceToParent;
		tmp = tmp->parent;
	}
	return distance;
}


//TODO - this is no longer in use really.
void PhylogeneticTree::fromNewick(const string& newick)
{
	//parse Newick!
	//start with a bracket, end with a bracket and a semicolon
	//TODO - check if the newick string is well formed
	bool endReached = false;
	unsigned int i = 0;
	unsigned int ids = 0;
	unsigned int sequenceNo;
	string seqName;
	//FIXME - write a regex that accepts any string
	regex regDescDist("((S[0-9]+):([0-9\\.]+)).*");
	//regex regDescDist("((.?+):([0-9\\.]+)).*");
	regex regInterDist(":([0-9\\.]+).*");
	Node *tmpNode, *tmpParent, *tmpCurrent;
	stack<Node*> workNodes;

	DEBUG("K-mer newick tree : " << newick);
	//cout << newick << endl;

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
			DUMP("newick (");
			tmpNode = new Node(++ids);
			nodes.push_back(tmpNode);
			workNodes.push(tmpNode);
			i++;
		}
		else if(newick[i] == ',')
		{
			DUMP("newick ,");
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

			DUMP("newick )");
			tmpCurrent = workNodes.top();
			workNodes.pop();
			if (workNodes.size() > 0)
			{
				tmpParent = workNodes.top();
				tmpParent->setChild(tmpCurrent);
				tmpCurrent->setParent(tmpParent);
			}


			i++;
		}
		else if(newick[i] == ';')
		{
			DUMP("newick ;");
			endReached = true;
			i++;
		}
		else
		{
			smatch basematch;
			string sstr = newick.substr(i);

			if (regex_match(sstr, basematch, regDescDist))
			{


				seqName = basematch[2].str();
				sequenceNo = inputSeqs->getSequenceId(seqName);
				currentDistance =  stof(basematch[3].str());

				//cerr << "Newick name, id and distance : " << seqName << " " << sequenceNo << " " << currentDistance << endl;

				tmpCurrent = workNodes.top();
				tmpCurrent->setSequenceNumber(sequenceNo);
				tmpCurrent->setDistance(currentDistance);
				tmpCurrent->setLeaf();
				this->leafNodes[sequenceNo] = tmpCurrent;
				i += basematch[1].str().size();
			}
			else if (regex_match(sstr, basematch, regInterDist))
			{
				currentDistance =   stof(basematch[1].str());

				//cerr << "Newick distance match: " << currentDistance << endl;

				tmpCurrent = workNodes.top();
				tmpCurrent->setDistance(currentDistance);
				i += basematch[1].str().size() + 1;
			}
		}
	}

	DEBUG("Newick tree parsed");

}

double PhylogeneticTree::distanceById(unsigned int n1, unsigned int n2) {
	Node* nd1 = leafNodes[n1];
	Node* nd2 = leafNodes[n2];
	Node* mrca = mostRecentAncestor(nd1,nd2);

	double distance = 0;

	for (Node* tmp : {nd1,nd2}){
		while (tmp != mrca){
			distance += tmp->distanceToParent;
			tmp = tmp->parent;
		}
	}
	return distance;

}

double PhylogeneticTree::distanceByName(string& n1, string& n2) {
}

} /* namespace EBC */
