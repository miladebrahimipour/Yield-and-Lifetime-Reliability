#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <algorithm>
#include "gate.h"
using namespace std;

#pragma once
struct nodes
{
	int ID;
	string gate;
	int numberOfInputs;
	vector<int> input;
	int out;
};


void ReadFromFile(vector<int> &inputs, vector<int> &outputs, vector<nodes> &node)
{
	int size = 0, size1 = 0, size2 = 0, size3 = 0;
	int i, o, n = 0, j = 0;
	string STRING;
	basic_string <char> substr;
	basic_string <char>::size_type indexCh1, indexCh2, indexCh3;
	fstream infile;
	infile.open("C7552.txt");
	while (!infile.eof()) // To get you all the lines.
	{
		getline(infile, STRING); // Saves the line in STRING.
		if (STRING[0] == '#' || STRING == "")
			continue;
		else if ((STRING.find("=")) != string::npos)
		{
			n = 0;
			size3 = 0;
			size++;
			indexCh2 = STRING.find_first_of("=");
			indexCh3 = STRING.find_first_of("(");
			substr = STRING.substr(0, indexCh2 - 1);
			node.resize(size);
			node[size - 1].ID = node[size - 1].out = stoi(substr);
			substr = STRING.substr(indexCh2 + 2, ((indexCh3 - 1) - (indexCh2 + 2) + 1));
			node[size - 1].gate = substr;
			//substr.clear();
			for (i = indexCh3, j = 0; i < STRING.length(); i++, j++)
			{
				if (STRING[i] == '(' || STRING[i] == ' ')
					continue;
				else if (STRING[i] == ',')
				{
					n++;
					node[size - 1].input.resize(n);
					node[size - 1].input[n - 1] = stoi(substr);
					size3 = 0;
					substr.clear();
				}
				else if (STRING[i] == ')')
				{
					n++;
					node[size - 1].input.resize(n);
					node[size - 1].input[n - 1] = stoi(substr);
					size3 = 0;
					substr.clear();
					break;
				}
				else
				{
					size3++;
					substr.resize(size3);
					substr[size3 - 1] = STRING[i];
				}

			}
			node[size - 1].numberOfInputs = node[size - 1].input.size();
		}
		else
		{
			indexCh1 = STRING.find("INPUT(");
			if (indexCh1 != string::npos)
			{
				size1++;
				indexCh2 = STRING.find_first_of("(");
				inputs.resize(size1);
				substr = STRING.substr(indexCh2 + 1, STRING.length());
				substr.pop_back();
				inputs[size1 - 1] = stoi(substr);
			}
			indexCh3 = STRING.find("OUTPUT(");
			if (indexCh3 != string::npos)
			{
				size2++;
				indexCh2 = STRING.find_first_of("(");
				outputs.resize(size2);
				substr = STRING.substr(indexCh2 + 1, STRING.length());
				substr.pop_back();
				outputs[size2 - 1] = stoi(substr);
			}
		}


	}
	infile.close();
}


void topologicalSort(const vector<int> &inputs, const vector<nodes> &node, vector<gate> &topological, 
	                 double Vdd, double Vth, double inputsSignalProbability, double activityFactor, double inputTransition)
{
	int i, j, temp;
	//*******************************   circuit generator variables
	int z;
	vector<int> t;
	t.resize(0);
	//*******************************
	int cnt = 0;
	int size = 0;
	bool flag = false;
	queue<int> q;
	for (i = 0; i < inputs.size(); i++)
		q.push(inputs[i]);

	vector<int> inputsNumber;
	inputsNumber.resize(node.size());
	for (i = 0; i < node.size(); i++)
	{
		inputsNumber[i] = node[i].numberOfInputs;
	}
	while (!q.empty())
	{
		temp = q.front();
		q.pop();
		size++;
		topological.resize(size);
		topological[size - 1].setID(temp);
		for (i = 0; i < node.size(); i++)
		{
			//*****************************   circuit generator
			if (node[i].ID == temp)  
			{
				topological[size - 1].set(node[i].ID, node[i].gate, node[i].numberOfInputs, node[i].input, node[i].out);
				topological[size - 1].setInputCapacitance();
				flag = true;
				continue;
			}
			//*****************************
			for (j = 0; j < node[i].numberOfInputs; j++)
			{
				if (node[i].input[j] == temp)
				{
					inputsNumber[i]--;
					if (inputsNumber[i] == 0)
						q.push(node[i].out);
				}
			}
		}
		cnt++;
		//*************************   circuit generator
		if (flag == false) 
		{
			for (z = 0; z < inputs.size(); z++)
			{
				if (inputs[z] == temp)
				{
					topological[size - 1].set(inputs[z], "INPUT", 0, t, inputs[z]);
					topological[size - 1].SetCircuiteInputs(Vdd, Vth, inputsSignalProbability, activityFactor);
					flag = false;
				}
			}
		}
		//*************************
	}
}


void SparsGenerator(vector<gate> &topological, const vector<nodes> &node)
{
	int i, j, k, z;
	int size;
	//spars.resize(topological.size());
	for (i = 0; i < topological.size(); i++)
	{
		size = 0;
		for (j = 0; j < node.size(); j++)
		{
			for (k = 0; k < node[j].numberOfInputs; k++)
			{
				if (node[j].input[k] == topological[i].getID())
				{
					size++;
					topological[i].output.resize(size);
					for (z = 0; z < topological.size(); z++)
					{
						if (topological[z].getID() == node[j].ID)
							topological[i].output[size - 1] = z;
					}
				}
			}
		}
	}

	//*********************************************************** Index Input creating
	for (i = topological.size() - 1; i >= 0; i--)
	{
		if (topological[i].getInputs().size() > 0)
		{
			for (j = 0; j < topological[i].getInputs().size(); j++)
			{
				for (k = 0; k < topological.size(); k++)
				{
					if (topological[i].getInputs()[j] == topological[k].getID())
						topological[i].getIndex_Inputs()[j] = k;
				}
			}
		}
	}
}


void UpdateNodes(vector<gate> &topological, double Vdd, double Vth)
{
	int i, j, k;
	double C = 0;
	for (i = 0; i < topological.size(); i++)
	{
		C = 0;
		if (topological[i].getGateName() != "INPUT")
			topological[i].Initialization(Vdd, Vth);
		for (j = 0; j < topological[i].output.size(); j++)
		{
			for (k = 0; k < topological[topological[i].output[j]].getInputs().size(); k++)
			{
				if (topological[i].getID() == topological[topological[i].output[j]].getInputs()[k])
				{
					topological[topological[i].output[j]].setSignalProbability(topological[i].getSignalProbability(), k);
					C += topological[topological[i].output[j]].getInputCapacitance(k);
					break;
				}
			}
		}
		topological[i].setCL(C);
	}
}


void ReadFromFile1(vector<int> &inputs, vector<int> &outputs, vector<nodes> &node, string FileName)
{
	int size = 0, size1 = 0, size2 = 0, size3 = 0;
	int i, o, n = 0, j = 0;
	string STRING;
	basic_string <char> substr;
	basic_string <char>::size_type indexCh1, indexCh2, indexCh3;
	fstream infile;
	infile.open(FileName);
	while (!infile.eof()) // To get you all the lines.
	{
		getline(infile, STRING); // Saves the line in STRING.
		if (STRING[0] == '#' || STRING == "")
			continue;
		else if ((STRING.find("=")) != string::npos)
		{
			n = 0;
			size3 = 0;
			size++;
			indexCh2 = STRING.find_first_of("=");
			indexCh3 = STRING.find_first_of("(");
			substr = STRING.substr(0, indexCh2 - 1);
			node.resize(size);
			node[size - 1].ID = node[size - 1].out = stoi(substr);
			substr = STRING.substr(indexCh2 + 2, ((indexCh3 - 1) - (indexCh2 + 2) + 1));
			node[size - 1].gate = substr;
			//substr.clear();
			for (i = indexCh3, j = 0; i < STRING.length(); i++, j++)
			{
				if (STRING[i] == '(' || STRING[i] == ' ')
					continue;
				else if (STRING[i] == ',')
				{
					n++;
					node[size - 1].input.resize(n);
					node[size - 1].input[n - 1] = stoi(substr);
					size3 = 0;
					substr.clear();
				}
				else if (STRING[i] == ')')
				{
					n++;
					node[size - 1].input.resize(n);
					node[size - 1].input[n - 1] = stoi(substr);
					size3 = 0;
					substr.clear();
					break;
				}
				else
				{
					size3++;
					substr.resize(size3);
					substr[size3 - 1] = STRING[i];
				}

			}
			node[size - 1].numberOfInputs = node[size - 1].input.size();
		}
		else
		{
			indexCh1 = STRING.find("INPUT(");
			if (indexCh1 != string::npos)
			{
				size1++;
				indexCh2 = STRING.find_first_of("(");
				inputs.resize(size1);
				substr = STRING.substr(indexCh2 + 1, STRING.length());
				substr.pop_back();
				inputs[size1 - 1] = stoi(substr);
			}
			indexCh3 = STRING.find("OUTPUT(");
			if (indexCh3 != string::npos)
			{
				size2++;
				indexCh2 = STRING.find_first_of("(");
				outputs.resize(size2);
				substr = STRING.substr(indexCh2 + 1, STRING.length());
				substr.pop_back();
				outputs[size2 - 1] = stoi(substr);
			}
		}


	}
	infile.close();
}


double PowerConsumption(vector<gate> &topological)
{
	int i, j, k;
	double power = 0, init_power = 0;
	for (i = 0; i < topological.size(); i++)
	{
		if (topological[i].getGateName() != "INPUT")
		{
			for (j = 0; j < topological[i].output.size(); j++)
			{
				for (k = 0; k < topological[topological[i].output[j]].getInputs().size(); k++)
				{
					power += (topological[topological[i].output[j]].getWL_ratio() * topological[topological[i].output[j]].getInputCapacitance(k));
					init_power += (2 * topological[topological[i].output[j]].getInputCapacitance(k));
				}
			}
		}
	}
	return ((abs(power - init_power)) / (init_power)) * 100;
}