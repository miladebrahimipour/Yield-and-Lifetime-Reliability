#include<iostream>
#include <vector>
#include<string>
using namespace std;


struct Joint_PDF
{
	double Aging_mean;
	double Criticality_mean;
	double Aging_Variance;
	double Criticality_Variance;
	double Correlation_coefficient;
};


#pragma once
class gate
{
public:
	gate();
	virtual ~gate();
	void signalProbability();
	void activityFactor();
	double gamma(int t, double sp);// double gamma(int t, double sp);
	double dVth(int t, double sp);
	void Initialization(double, double);
	double maxSignalProbability();
	void SetCircuiteInputs(double vdd, double vth, double sp, double a);
	void setSignalProbability(double sp, int i);
	void setActivityFactor(double a, int i);
	double getSignalProbability();
	double getActivityFactor();
	vector<int>& getInputs();
	vector<int>& getIndex_Inputs();
	void setInputCapacitance();
	double getInputCapacitance(int i);
	void setCL(double c);
	void set(int I, string g, int numberOfInputs, vector<int> input, int out);
	void setID(int I);
	int getID();
	double getCL();
	string getGateName();
	void setOutputProbability(double i);
	void setOutSwitchingActivity(double i);
	//***********************************************************
	// mean
	void setDelay_mean();
	double getDelay_mean();

	void setAging_Delay_mean(int t);
	double getAging_Delay_mean();

	void setArrival_mean(double i);
	double getArrival_mean();

	void setRAT_mean(double i);
	double getRAT_mean();

	void setCriticality_mean(double i);
	double getCriticality_mean();

	// Vth sensitivity
	void setDelay_Vth_sensitivity(int i);
    double getDelay_Vth_sensitivity(int i);

	void setAging_Delay_Vth_sensitivity(int i, int t);//void setAging_Delay_Vth_sensitivity(int i, int t);
	double getAging_Delay_Vth_sensitivity(int i);

	void setArrival_Vth_sensitivity(int i, double j);
	double getArrival_Vth_sensitivity(int i);

	void setRAT_Vth_sensitivity(int i, double j);
	double getRAT_Vth_sensitivity(int i);

	void setCriticality_Vth_sensitivity(int i, double j);
	double getCriticality_Vth_sensitivity(int i);
	vector<double>& getCriticality_Vth_sensitivity();

	// Vth rand sensitivity
	void setDelay_Vth_rand_sensitivity();
	double getDelay_Vth_rand_sensitivity();

	void setAging_Delay_Vth_rand_sensitivity(int t);//void setAging_Delay_Vth_rand_sensitivity(int t);
	double getAging_Delay_Vth_rand_sensitivity();

	void setArrival_Vth_rand_sensitivity(double i);
	double getArrival_Vth_rand_sensitivity();

	void setRAT_Vth_rand_sensitivity(double i);
	double getRAT_Vth_rand_sensitivity();

	void setCriticality_Vth_rand_sensitivity(double i);
	double getCriticality_Vth_rand_sensitivity();

	//***************************************************************
	vector<double>& getArrival_Vth_sensitivity();
	vector<double>& getRAT_Vth_sensitivity();
	vector<double>& getAging_Delay_Vth_sensitivity();

	void setCriticality(double i);
	double getCriticality();

	void setAverage_Criticality(double i);
	double getAverage_Criticality();

	void setAging(double i);
	double getAging();

	void setRank(double i);
	double getRank();

	void setWL_ratio(double i);
	double getWL_ratio();

	void setJoint(Joint_PDF i);
	Joint_PDF getJoint_PDF();

	void setCritical_Aging(double i);
	double getCritical_Aging();

	vector<int> output;

private:
	double Vth;
	double Vdd;
	int ID;
	string gateName;
	int numberOfInputs;
	vector<int> input;
	vector<int> index_input;    //**********************************************************************************************************  Index_Input
	int out;
	vector<double> inputProbability;
	double outputProbability;
	double switchingActivity;
	vector<double> Cin;
	double CL;
	//double dVth;

	double Delay_mean;
	double Aging_Delay_mean;
	double Arrival_mean;
	double RAT_mean;
	double Criticality_mean;

	vector<double> Delay_Vth_sensitivity;
	vector<double> Aging_Delay_Vth_sensitivity;
	vector<double> Arrival_Vth_sensitivity;
	vector<double> RAT_Vth_sensitivity;
	vector<double> Criticality_Vth_sensitivity;

	double Delay_Vth_rand_sensitivity;
	double Aging_Delay_Vth_rand_sensitivity;
	double Arrival_Vth_rand_sensitivity;
	double RAT_Vth_rand_sensitivity;
	double Criticality_Vth_rand_sensitivity;

	double Criticality;
	double Average_Criticality;
	double Aging;
	double Rank;
	double WL_ratio;
	Joint_PDF Joint;
	double Critical_Aging;
};