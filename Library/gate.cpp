#include "stdafx.h"
#include "gate.h"
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include "PreDefinedTables.h"
using namespace std;


gate::gate()
{
}


gate::~gate()
{
}


void gate::signalProbability()
{
	double temp1 = 0;
	outputProbability = 1;
	if (gateName == "AND")
	{
		for (int i = 0; i < numberOfInputs; i++)
			outputProbability *= inputProbability[i];
	}
	else if (gateName == "OR")
	{
		for (int i = 0; i < numberOfInputs; i++)
		{
			temp1 *= (1 - inputProbability[i]);
		}
		outputProbability = 1 - temp1;
	}
	else if (gateName == "XOR")
	{
		outputProbability = ((inputProbability[0] * (1 - inputProbability[1])) + ((1 - inputProbability[0]) * inputProbability[1]));
	}
	else if (gateName == "NAND")
	{
		for (int i = 0; i < numberOfInputs; i++)
			outputProbability *= inputProbability[i];
		outputProbability = 1 - outputProbability;
	}
	else if (gateName == "NOR")
	{
		for (int i = 0; i < numberOfInputs; i++)
		{
			temp1 *= (1 - inputProbability[i]);
		}
		outputProbability = temp1;
	}
	else if (gateName == "NOT")
		outputProbability = 1 - inputProbability[0];
	else if (gateName == "BUFF")
		outputProbability = inputProbability[0];
	else
		return;
}


void gate::activityFactor()
{
	switchingActivity = outputProbability * (1 - outputProbability);
}


double gate::gamma(int t , double sp)
{
	double temp, temp1, temp2, temp3;
	const double q = 1.6e-19;
	const double Tox = 1.5e-9;
	const double T = 300;
	const double Ea = 0.49;
	const double To = 1e10;
	const double Eo = 0.335e9;
	const double K = 9.4868e25;
	const double k = 8.617385E-5;
	const double epsox = 3.9 * 8.8542e-12;
	const double e1 = 0.9e-9;
	const double e2 = 0.5e-4;
	const double te = 1.82e-9;//Tox; 
	const double n = 1.0 / 6;
	double Cox = epsox / Tox;
	double Eox = (Vdd - Vth) / Tox;
	long double C = pow(To, (-1)) * exp(-(Ea / (k * T)));
	double Kv = pow(((q * Tox) / epsox), 3) * pow(K, 2) * Cox * (Vdd - Vth) * sqrt(C) * exp((2 * Eox) / Eo);
	/*if (sp == 1)
		temp3 = 2 * n * pow(Kv, ((2 * n) - 1)) * pow(t, n) * (Kv / (Vdd - Vth));
	else
	{
		long double Bt = 1 - (((2 * e1 * te) + sqrt(e2 * C * (1 - sp) * Tclk)) / ((2 * Tox) + sqrt(C * t)));
		//temp = (pow(Kv, 2) * sp * Tclk);
		temp = 2 * n * pow(Kv, ((2 * n) - 1)) * pow(sp, n) * pow(t, n);
		temp1 = pow((1 - pow(Bt, 1 / (2 * n))), (2 * n));
		temp2 = temp / temp1;
		temp3 = temp2 * (Kv / (Vdd - Vth));//2 * n * pow(Kv, ((2 * n) - 1)) * pow(t, n) * pow(sp, n)* (Kv / (Vdd - Vth));
		//dVth = pow((sqrt(temp) / temp2), 2 * n);
	}*/
	if (sp == 1)
	{
		temp3 = 2 * n * pow(Kv, ((2 * n) - 1)) * pow(t , n) * (Kv / (Vdd - Vth));
	//	dVth = pow((temp * t), n);
	}
	else
	{
		temp = pow((pow(n, 2) * sp * C * t * 0.001), n);
		temp1 = temp * 2 * n * pow(Kv, ((2 * n) - 1));
		temp2 = pow((0.81 * pow(Tox, 2) * (1 - sp)), n);
		temp3 = (temp1 / temp2) * (Kv / (Vdd - Vth));
	}
	return temp3;
	//*********************************************************************
	/*
	if (sp == 1)
	{
		temp = pow(Kv, 2);
		dVth = pow((temp * t), n);
	}

	else
	{
		//long double Bt = 1 - (((2 * e1 * te) + sqrt(e2 * C * (1 - sp) * Tclk)) / ((2 * Tox) + sqrt(C * t)));
		//********************************************************************************************************
		long double temp = pow(Kv, 2) * pow(n, 2) * sp * C * t * 0.001;
		long double temp2 = 0.81 * pow(Tox, 2) * (1 - sp);
		dVth = pow(temp / temp2, n);
		//double g = (sp > (1 - sp)) ? 1 - sp : sp;
		//temp = (pow(Kv, 2) * sp * Tclk);// / g;
		//temp2 = 1 - pow(Bt, 1 / (2 * n));
		//dVth = pow((sqrt(temp) / temp2), 2 * n);
	}
	*/
}


double gate::dVth(int t, double sp)
{
	double temp, temp1, temp2, temp3;
	const double q = 1.6e-19;
	const double Tox = 1.5e-9;
	const double T = 300;
	const double Ea = 0.49;
	const double To = 1e10;
	const double Eo = 0.335e9;
	const double K = 9.4868e25;
	const double k = 8.617385E-5;
	const double epsox = 3.9 * 8.8542e-12;
	const double e1 = 0.9e-9;
	const double e2 = 0.5e-4;
	const double te = 1.82e-9;//Tox; 
	const double n = 1.0 / 6;
	double Cox = epsox / Tox;
	double Eox = (Vdd - Vth) / Tox;
	long double C = pow(To, (-1)) * exp(-(Ea / (k * T)));
	double Kv = pow(((q * Tox) / epsox), 3) * pow(K, 2) * Cox * (Vdd - Vth) * sqrt(C) * exp((2 * Eox) / Eo);
	/*if (sp == 1)
	temp3 = 2 * n * pow(Kv, ((2 * n) - 1)) * pow(t, n) * (Kv / (Vdd - Vth));
	else
	{
	long double Bt = 1 - (((2 * e1 * te) + sqrt(e2 * C * (1 - sp) * Tclk)) / ((2 * Tox) + sqrt(C * t)));
	//temp = (pow(Kv, 2) * sp * Tclk);
	temp = 2 * n * pow(Kv, ((2 * n) - 1)) * pow(sp, n) * pow(t, n);
	temp1 = pow((1 - pow(Bt, 1 / (2 * n))), (2 * n));
	temp2 = temp / temp1;
	temp3 = temp2 * (Kv / (Vdd - Vth));//2 * n * pow(Kv, ((2 * n) - 1)) * pow(t, n) * pow(sp, n)* (Kv / (Vdd - Vth));
	//dVth = pow((sqrt(temp) / temp2), 2 * n);
	}*/
	if (sp == 1)
	{
		temp = pow(Kv, 2);
		temp3 = pow((temp * t), n);
		//	dVth = pow((temp * t), n);
	}
	else
	{
		temp = pow(Kv, 2) * pow(n, 2) * sp * C * t * 0.001;
	    temp2 = 0.81 * pow(Tox, 2) * (1 - sp);
		//temp2 = pow((0.81 * pow(Tox, 2) * (1 - sp)), n);
		temp3 = pow(temp / temp2, n);
	}
	return temp3;
	//*********************************************************************
	/*
	if (sp == 1)
	{
	temp = pow(Kv, 2);
	dVth = pow((temp * t), n);
	}

	else
	{
	//long double Bt = 1 - (((2 * e1 * te) + sqrt(e2 * C * (1 - sp) * Tclk)) / ((2 * Tox) + sqrt(C * t)));
	//********************************************************************************************************
	long double temp = pow(Kv, 2) * pow(n, 2) * sp * C * t * 0.001;
	long double temp2 = 0.81 * pow(Tox, 2) * (1 - sp);
	dVth = pow(temp / temp2, n);
	//double g = (sp > (1 - sp)) ? 1 - sp : sp;
	//temp = (pow(Kv, 2) * sp * Tclk);// / g;
	//temp2 = 1 - pow(Bt, 1 / (2 * n));
	//dVth = pow((sqrt(temp) / temp2), 2 * n);
	}
	*/
}


void gate::Initialization(double vdd, double vth)
{
	Vdd = vdd;
	Vth = vth;
	signalProbability();
	activityFactor();
}


double gate::maxSignalProbability()
{
	double sp;
	sp = 1 - *min(inputProbability.begin(), inputProbability.end());
	return sp;
}


void gate::SetCircuiteInputs(double vdd, double vth, double sp, double a)
{
	Vdd = vdd;
	Vth = vth;
	outputProbability = sp;
	switchingActivity = a;
	Delay_mean = 0;
	Delay_Vth_rand_sensitivity = 0;
	Arrival_mean = 0;
	Arrival_Vth_rand_sensitivity = 0;
}

//**********************************************************************************************************  Index_Input
void gate::set(int I, string g, int n, vector<int> in, int o)   
{
	ID = I;
	gateName = g;
	numberOfInputs = n;
	Criticality = 0;
	WL_ratio = 2;
	if (n != 0)
	{
		input.resize(n);
		index_input.resize(n);         //************************************* Index_Input
		inputProbability.resize(n);
		Cin.resize(numberOfInputs);
		input = in;
	}
	if (gateName != "INPUT")
	{
		Delay_Vth_sensitivity.resize(16);
		Aging_Delay_Vth_sensitivity.resize(16);
		Arrival_Vth_sensitivity.resize(16);
		RAT_Vth_sensitivity.resize(16);
		Criticality_Vth_sensitivity.resize(16);
	}
	out = o;
}


void gate::setID(int I)
{
	ID = I;
}


int gate::getID()
{
	return ID;
}


void gate::setSignalProbability(double sp,int i)
{
	inputProbability[i] = sp;
}


double gate::getSignalProbability()
{
	return outputProbability;
}


double gate::getActivityFactor()
{
	return switchingActivity;
}


vector<int>& gate::getInputs()
{
	return input;
}

//**********************************************************************************************************  Index_Input
vector<int>& gate::getIndex_Inputs()
{
	return index_input;
}


void gate::setInputCapacitance()
{
	if (gateName == "AND")
	{
		switch (numberOfInputs)
		{
		case 2:
			Cin[0] = 0.9181;
			Cin[1] = 0.9746;
			break;
		case 3:
			Cin[0] = 0.8797;
			Cin[1] = 0.9275;
			Cin[2] = 0.9648;
			break;
		//case 4:
		//	Cin[0] = 0.8565;
		//	Cin[1] = 0.9023;
		//	Cin[2] = 0.9241;
		//	Cin[3] = 0.9445;
		//	break;
		default:
			Cin[0] = 0.8565;
			Cin[1] = 0.9023;
			Cin[2] = 0.9241;
			Cin[3] = 0.9445;
			for (int i = 4; i < numberOfInputs; i++)
				Cin[i] = 0.9445;
		}
	}

	else if (gateName == "OR")
	{
		switch (numberOfInputs)
		{
		case 2:
			Cin[0] = 0.9468;
			Cin[1] = 0.9419;
			break;
		case 3:
			Cin[0] = 0.9591;
			Cin[1] = 0.9401;
			Cin[2] = 0.9216;
			break;
		//case 4:
		//	Cin[0] = 0.9412;
		//	Cin[1] = 0.9383;
		//	Cin[2] = 0.9238;
		//	Cin[3] = 0.9142;
		//	break;
		default:
			Cin[0] = 0.9412;
			Cin[1] = 0.9383;
			Cin[2] = 0.9238;
			Cin[3] = 0.9142;
			for (int i = 4; i < numberOfInputs; i++)
				Cin[i] = 0.9142;
		}
	}

	else if (gateName == "NAND")
	{
		switch (numberOfInputs)
		{
		case 2:
			Cin[0] = 1.5990;
			Cin[1] = 1.6642;
			break;
		case 3:
			Cin[0] = 1.5903;
			Cin[1] = 1.6212;
			Cin[2] = 1.6504;
			break;
		//case 4:
			//Cin[0] = 1.5221;
			//Cin[1] = 1.5952;
			//Cin[2] = 1.6381;
			//Cin[3] = 1.6599;
		//	break;
		default:
			Cin[0] = 1.5221;
			Cin[1] = 1.5952;
			Cin[2] = 1.6381;
			Cin[3] = 1.6599;
			for (int i = 4; i < numberOfInputs; i++)
				Cin[i] = 1.6599;
		}
	}

	else if (gateName == "NOR")
	{
		switch (numberOfInputs)
		{
		case 2:
			Cin[0] = 1.7145;
			Cin[1] = 1.6513;
			break;
		case 3:
			Cin[0] = 1.7636;
			Cin[1] = 1.6638;
			Cin[2] = 1.6163;
			break;
			//case 4:
			//Cin[0] = 1.7368;
			//Cin[1] = 1.6741;
			//Cin[2] = 1.6360;
			//Cin[3] = 1.6059;
			//	break;
		default:
			Cin[0] = 1.7368;
			Cin[1] = 1.6741;
			Cin[2] = 1.6360;
			Cin[3] = 1.6059;
			for (int i = 4; i < numberOfInputs; i++)
				Cin[i] = 1.6059;
		}
	}

	else if (gateName == "NOT")
	{
		Cin[0] = 1.7002;
	}

	else if (gateName == "BUFF")
	{
		Cin[0] = 0.9747;
	}

	else if (gateName == "XOR")
	{
		Cin[0] = 2.2321;
		Cin[1] = 2.4115;
	}
}


double gate::getInputCapacitance(int i)
{
	return Cin[i];
}


void gate::setCL(double c)
{
	CL = c;
}


double gate::getCL()
{
	return CL;
}


string gate::getGateName()
{
	return gateName;
}


void gate::setOutputProbability(double i)
{
	outputProbability = i;
}


void gate::setOutSwitchingActivity(double i)
{
	switchingActivity = i;
}


//********************************************************* set and get functions for statistical parameters

// mean
void gate::setDelay_mean()
{
	double alpha = 1.5;// 0.1323;
	double Tox = 9.2e-010;
	double epsox = 3.9 * 8.8542e-12;
	double Cox = epsox / Tox;
	double mu = 0.00391;
	double beta = mu * Cox * WL_ratio;
	double temp = CL * Vdd;
	double temp1 = pow((Vdd - Vth), alpha);
	double temp2 = beta * temp1;
	Delay_mean = temp / temp2;
	Delay_mean *= pow(10, -15);
}


double gate::getDelay_mean()
{
	return Delay_mean;
}


void gate::setAging_Delay_mean(int t)
{
	double A = 8.08e-3;//0.003606994884410881;// 0.003214779955;// 0.002970042;// 8.08e-6;// 0.7648e-3; // 45.58e-3;//29.08e-3;
	double n = 1.0 / 6;
	double alpha = 1.5;
	double temp = alpha * Delay_mean;
	double temp1 = Vdd - Vth;
	double beta = temp / temp1;
	double temp2 = A * pow(maxSignalProbability(), n) * pow(t, n);
	Aging_Delay_mean = A * pow(maxSignalProbability(), n) * pow(t, n) * beta;// * dVth(10 * 365 * 24 * 60 * 60, maxSignalProbability());
}


double gate::getAging_Delay_mean()
{
	return Aging_Delay_mean;
}


void gate::setArrival_mean(double i)
{
	Arrival_mean = i;
}


double gate::getArrival_mean()
{
	return Arrival_mean;
}


void gate::setRAT_mean(double i)
{
	RAT_mean = i;
}


double gate::getRAT_mean()
{
	return RAT_mean;
}


void gate::setCriticality_mean(double i)
{
	Criticality_mean = i;
}


double gate::getCriticality_mean()
{
	return Criticality_mean;
}


// Vth sensitivity
void gate::setDelay_Vth_sensitivity(int i)
{
	double alpha = 1.5;// 0.1323;
	double temp = alpha * Delay_mean;
	double temp1 = Vdd - Vth;
	Delay_Vth_sensitivity[i] = temp / temp1;
}


double gate::getDelay_Vth_sensitivity(int i)
{
	return Delay_Vth_sensitivity[i];
}


void gate::setAging_Delay_Vth_sensitivity(int i, int t)
{
	double A = 8.08e-3; //0.003606994884410881;// 0.003214779955;// 0.002970042;// 0.003214779955;// 0.7648e-3;// 45.58e-3;
	double n = 1.0 / 6;
	double alpha = 1.5;// 0.1323;
	double temp = alpha * Delay_mean;
	double temp1 = Vdd - Vth;
	double beta = temp / temp1;
	double gama = 0.1423250;//-(0.14532);// gamma(10 * 365 * 24 * 60 * 60, maxSignalProbability());
	Aging_Delay_Vth_sensitivity[i] = A * pow(maxSignalProbability(), n) * pow(t, n) * beta * gama;
}


double gate::getAging_Delay_Vth_sensitivity(int i)
{
	return Aging_Delay_Vth_sensitivity[i];
}


void gate::setArrival_Vth_sensitivity(int i, double j)
{
	Arrival_Vth_sensitivity[i] = j;
}


double gate::getArrival_Vth_sensitivity(int i)
{
	return Arrival_Vth_sensitivity[i];
}


void gate::setRAT_Vth_sensitivity(int i, double j)
{
	RAT_Vth_sensitivity[i] = j;
}


double gate::getRAT_Vth_sensitivity(int i)
{
	return RAT_Vth_sensitivity[i];
}


void gate::setCriticality_Vth_sensitivity(int i, double j)
{
	Criticality_Vth_sensitivity[i] = j;
}


double gate::getCriticality_Vth_sensitivity(int i)
{
	return Criticality_Vth_sensitivity[i];
}


vector<double>& gate::getCriticality_Vth_sensitivity()
{
	return Criticality_Vth_sensitivity;
}


// Vth rand sensitivity
void gate::setDelay_Vth_rand_sensitivity()
{
	double alpha = 1.5;// 0.1323;
	double temp = alpha * Delay_mean;
	double temp1 = Vdd - Vth;
	Delay_Vth_rand_sensitivity = temp / temp1;
}


double gate::getDelay_Vth_rand_sensitivity()
{
	return Delay_Vth_rand_sensitivity;
}


void gate::setAging_Delay_Vth_rand_sensitivity(int t)
{
	double A = 8.08e-3;//0.003606994884410881;// 0.003214779955;// 0.002970042;// 0.0032147799557;// 0.7646e-3;
	double n = 1.0 / 6;
	double alpha = 1.5;// 0.1323;
	double temp = alpha * Delay_mean;
	double temp1 = Vdd - Vth;
	double beta = temp / temp1;
	double gama = 0.1423250;//-(0.14532);// gamma(10 * 365 * 24 * 60 * 60, maxSignalProbability());
	Aging_Delay_Vth_rand_sensitivity = A * pow(maxSignalProbability(), n) * pow(t, n) * gama * beta;
}


double gate::getAging_Delay_Vth_rand_sensitivity()
{
	return Aging_Delay_Vth_rand_sensitivity;
}


void gate::setArrival_Vth_rand_sensitivity(double i)
{
	Arrival_Vth_rand_sensitivity = i;
}


double gate::getArrival_Vth_rand_sensitivity()
{
	return Arrival_Vth_rand_sensitivity;
}


void gate::setRAT_Vth_rand_sensitivity(double i)
{
	RAT_Vth_rand_sensitivity = i;
}


double gate::getRAT_Vth_rand_sensitivity()
{
	return RAT_Vth_rand_sensitivity;
}


void gate::setCriticality_Vth_rand_sensitivity(double i)
{
	Criticality_Vth_rand_sensitivity = i;
}


double gate::getCriticality_Vth_rand_sensitivity()
{
	return Criticality_Vth_rand_sensitivity;
}


//******************************************************************
vector<double>& gate::getArrival_Vth_sensitivity()
{
	return Arrival_Vth_sensitivity;
}


vector<double>& gate::getRAT_Vth_sensitivity()
{
	return RAT_Vth_sensitivity;
}


vector<double>& gate::getAging_Delay_Vth_sensitivity()
{
	return Aging_Delay_Vth_sensitivity;
}


void gate::setCriticality(double i)
{
	Criticality = i;
}


double gate::getCriticality()
{
	return Criticality;
}


void gate::setAverage_Criticality(double i)
{
	Average_Criticality = i;
}


double gate::getAverage_Criticality()
{
	return Average_Criticality;
}


void gate::setAging(double i)
{
	Aging = i;
}


double gate::getAging()
{
	return Aging;
}


void gate::setRank(double i)
{
	Rank = i;
}


double gate::getRank()
{
	return Rank;
}


void gate::setWL_ratio(double i)
{
	WL_ratio = i;
}


double gate::getWL_ratio()
{
	return WL_ratio;
}


void gate::setJoint(Joint_PDF i)
{
	Joint.Aging_mean = i.Aging_mean;
	Joint.Criticality_mean = i.Criticality_mean;
	Joint.Aging_Variance = i.Aging_Variance;
	Joint.Criticality_Variance = i.Criticality_Variance;
	Joint.Correlation_coefficient = i.Correlation_coefficient;
}


Joint_PDF gate::getJoint_PDF()
{
	return Joint;
}


void gate::setCritical_Aging(double i)
{
	Critical_Aging = i;
}


double gate::getCritical_Aging()
{
	return Critical_Aging;
}
