// SA_Optimization.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <ctime>
#include <iomanip>
#include "gate.h"
#include "Functions.h"
#include "Timing_Analysis.h"
#include "Optimization.h"
using namespace std;


const double Vdd = 1.1;
const double Vth = 0.4;
const double inputsSignalProbability = 0.5;
const double activityFactor = 0.5;
const double inputTransition = 0.067257142857;


void Initialization1(vector<int> &, vector<nodes> &, vector<int> &, vector<gate> &, string);
int _tmain(int argc, _TCHAR* argv[])
{
	int counter = 600000;
	int i, j, count = 0, count1;
	double temp = 0, temp1;
	double Constraint;
	double Area = 0, Area1 = 0, delta_Area, delta_Area1;
	double mean, mean1, mean_sub, var, var1, var_sub, yield, lifetime_yield;
	vector<int> inputs, outputs;
	vector<nodes> node;
	vector<PartialNode> PartialNodes;
	vector<gate> topological;
	vector<Critical> Criticality, Criticality1;
	vector<Age_Critical> Aging_Criticality;
	string FileName, output_FileName, outfile1;
	Sink sink;
	fstream Results, Results1;
	clock_t clk;
	srand(time(0));
	output_FileName = "Revision_Power_Without_Yield.txt";
	outfile1 = "worst_case_delay_ISCAS.txt";
	FileName = "multiplier.bench";

	Results.open(output_FileName, ios::out | ios::ate | ios::app);
	Results1.open(outfile1, ios::out | ios::app);
	clk = clock();
	Initialization1(inputs, node, outputs, topological, FileName);
	for (i = 0; i < topological.size(); i++)
	{
		if (topological[i].output.size() == 0)
			topological[i].setCL(10);
	}
	Spatial_Correlation(topological, inputs);
	SSTA(topological, &sink, 0.6, 1.1);
	//*****************************************************************************
	temp = Variance(sink.Arrival_Vth_sensitivity, sink.Arrival_Vth_rand_sensitivity);
	Results << left << setw(21) << "Bench" << setw(21) << "initial_mean" << setw(21) << "initial_variance" << setw(21) << "lifetime_mean" << setw(21)
		                              << "lifetime_variance" << setw(21) << "opt_lifetime_mean" << setw(21) << "opt_lifetime_variance" << setw(21) << "     Area" << endl;
	Results << left << setw(21) << FileName << setw(21) << sink.Arrival_mean << setw(21) << temp;
	cout << sink.Arrival_mean << "\t" << temp << endl;
	Results1 << "Bench: " << FileName << endl << "initial_delay: " << sink.Arrival_mean + (3 * sqrt(temp)) << endl;
	//*****************************************************************************
	Criticality_Computation(topological, &sink);
	Rank_Criticality(topological, Criticality);
	Constraint = temp1 = Yield_Computation(sink.Arrival_mean, sink.Arrival_Vth_sensitivity, sink.Arrival_Vth_rand_sensitivity, sink.RequiredArrivalTime);
	//Yield_Optimization(topological, Criticality, &sink, &count, 0.95, temp1);
	mean = sink.Arrival_mean;
	var = Variance(sink.Arrival_Vth_sensitivity, sink.Arrival_Vth_rand_sensitivity);
	yield = Yield_Computation(sink.Arrival_mean, sink.Arrival_Vth_sensitivity, sink.Arrival_Vth_rand_sensitivity, sink.RequiredArrivalTime);;
	/*for (i = 0; i < topological.size(); i++)
	{
		if (topological[i].getGateName() != "INPUT")
		{
			Area1 += topological[i].getWL_ratio();
			Area += 2;
		}
	}
	delta_Area = ((Area1 * 100) / Area) - 100;*/
	SSTA(topological, &sink, 0.6, 1.1, lifetime);
	//*****************************************************************************
	temp = 0;
	temp = Variance(sink.Arrival_Vth_sensitivity, sink.Arrival_Vth_rand_sensitivity);
	Results << left << setw(21) << sink.Arrival_mean << setw(21) << temp;
	Results1 << "lifetime_delay: " << sink.Arrival_mean + (3 * sqrt(temp)) << endl;
	cout << "**************************************************************************************" << endl;
	cout << sink.Arrival_mean << "\t" << temp << endl;
	//*****************************************************************************
	//Criticality_Computation(topological, &sink);
	//Aging_Computation(topological, &sink);
	//Pareto_Ranking(topological, &sink, Aging_Criticality);
	//Pareto_Ranking_Multiply(topological, &sink, Criticality1);
	Critical_Aging_Computation(topological, &sink);
	Rank_Critical_Aging(topological, Criticality1);
	//Constraint = 1.05 * (sink.Arrival_mean + (3 * sqrt(Variance(sink.Arrival_Vth_sensitivity, sink.Arrival_Vth_rand_sensitivity))));//sink.delay_mean + (3 * sqrt(Variance(sink.delay_Vth_sensitivity, sink.delay_Vth_rand_sensitivity)));
	//Constraint = 1.1 * sink.RequiredArrivalTime;
	Constraint = sink.Lifetime_RequiredArrivalTime;
	temp1 = sink.Arrival_mean + (3 * sqrt(Variance(sink.Arrival_Vth_sensitivity, sink.Arrival_Vth_rand_sensitivity)));
	//GB_Gate_Sizing(topological, Aging_Criticality, inputs, &sink, &count1, Constraint, temp1);
	//GB_Gate_Sizing_Multiply(topological, Criticality1, inputs, &sink, &count1, Constraint, temp1);
	GB_Gate_Sizing_Critical_Aging(topological, Criticality1, inputs, &sink, &count1, Constraint, temp1);
	//GB_Gate_Sizing_Critical_Aging_itr(topological, Criticality1, inputs, &sink, &count1, Constraint, temp1);
	lifetime_yield = Yield_Computation(sink.Arrival_mean, sink.Arrival_Vth_sensitivity, sink.Arrival_Vth_rand_sensitivity, Constraint);
	mean1 = sink.Arrival_mean;
	var1 = Variance(sink.Arrival_Vth_sensitivity, sink.Arrival_Vth_rand_sensitivity);
	clk = clock() - clk;
	for (i = 0; i < topological.size(); i++)
	{
		if (topological[i].getGateName() != "INPUT")
		{
			Area1 += topological[i].getWL_ratio();
			Area += 2;
		}
	}
	delta_Area1 = ((Area1 - Area) / Area) * 100;
	//****************************************************************************
	temp = 0;
	temp = Variance(sink.Arrival_Vth_sensitivity, sink.Arrival_Vth_rand_sensitivity);
	Results << left << setw(21) << sink.Arrival_mean << setw(21) << temp << setw(21) << "     " << "Area : " << delta_Area1 << endl;
	Results << "Runtime : " << ((double)clk) / (CLOCKS_PER_SEC) << endl;
	Results << "Yield_Count : " << count << "                     " << "Lifetime_Count : " << count1 << endl;
	Results1 << "opt_lifetime_delay: " << sink.Arrival_mean + (3 * sqrt(temp)) << endl;
	//****************************************************************************
	SSTA(topological, &sink, 0.6, 1.1, initial_optimization);
	temp = Variance(sink.Arrival_Vth_sensitivity, sink.Arrival_Vth_rand_sensitivity);
	temp1 = Yield_Computation(sink.Arrival_mean, sink.Arrival_Vth_sensitivity, sink.Arrival_Vth_rand_sensitivity, sink.RequiredArrivalTime);
	SUB_Delay_Distribution(mean1, var1, sink.Arrival_mean, temp, mean_sub, var_sub);
	Results << "Final_Opt_Circuit_mean: " << sink.Arrival_mean << endl << "Final_Opt_Circuit_Var: " << temp << endl << "Final_Opt_Circuit_Yield: " << temp1 << endl;
	Results << "Yield_Opt_mean: " << mean << endl << "Yield_Opt_var: " << var << endl << "Yield_yield: " << yield << endl << "GB_yield: " << lifetime_yield << endl;
	Results << "RequiredArrivalTime: " << sink.RequiredArrivalTime << "Lifetime_Required_Arrival_Time: " << Constraint << endl;
	Results << "GB_mean: " << mean_sub << endl << "GB_var: " << var_sub << endl;
	Results << "Power: " << PowerConsumption(topological) << endl;
	//****************************************************************************
	Results << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	Results1 << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	//****************************************************************************
	Results.close();
	return 0;
}


void Initialization1(vector<int> &inputs, vector<nodes> &node, vector<int> &outputs, vector<gate> &topological, string FileName)
{
	ReadFromFile1(inputs, outputs, node, FileName);
	topologicalSort(inputs, node, topological, Vdd, Vth, inputsSignalProbability, activityFactor, inputTransition);
	SparsGenerator(topological, node);
	UpdateNodes(topological, Vdd, Vth);
}