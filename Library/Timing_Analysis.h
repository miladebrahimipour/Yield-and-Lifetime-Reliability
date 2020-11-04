#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <Windows.h>
#include <cstdlib>
#include "Functions.h"
#include "gate.h"
#include "engine.h"
using namespace std;


#pragma once
struct Critical
{
	double criticality;
	int index;
};


struct Age_Critical
{
	int Pareto_index;
	int index;
};


enum mode { initial = 0, initial_optimization, lifetime, lifetime_optimization };


double Normal_Distribution_CDF(double x)
{
	double temp;
	int temp1, i, j = 0;
	bool flag = false;
	string STRING, str_temp, s, tz;
	basic_string <char>::size_type indexCh1, indexCh2;
	fstream infile;
	infile.open("Normal_Distribution_Predefined_table.txt");
	str_temp = to_string(x);
	str_temp = str_temp.substr(0, 4);
	temp = stod(str_temp);
	str_temp.erase(3, str_temp.length());
	if (x >= (temp + 0.005))
		temp += 0.01;
	s = to_string(temp); 
	tz = s.substr(3, 1);
	j = stoi(to_string(temp).substr(3, 1));
	while (!infile.eof()) // To get you all the lines.
	{
		getline(infile, STRING); // Saves the line in STRING.
		if (STRING.substr(0, 3) == str_temp)
		{
			STRING.erase(0, 3);
			i = 0;
			while (true)
			{
				indexCh1 = STRING.find_first_of(",");
				STRING.erase(indexCh1, 1);
				indexCh2 = STRING.find_first_of(",");
				if (indexCh2 == string::npos)
					indexCh2 = STRING.length();
				if (i == j)
					return stod(STRING.substr(indexCh1, indexCh2 - indexCh1));
				STRING.erase(0, indexCh2);
				i++;
			}
		}
	}
	infile.close();
}


void swap(vector<Critical> &Criticality, int i, int j)
{
	double temp;
	int temp1;

	temp = Criticality[i].criticality;
	temp1 = Criticality[i].index;
	Criticality[i].criticality = Criticality[j].criticality;
	Criticality[i].index = Criticality[j].index;
	Criticality[j].criticality = temp;
	Criticality[j].index = temp1;
}


double Phi(double x)
{
	if (x <= 0.05)
		return 0.01994;
	else if (x <= 0.10)
		return 0.03983;
	else if (x <= 0.15)
		return 0.05962;
	else if (x <= 0.20)
		return 0.07926;
	else if (x <= 0.25)
		return 0.09871;
	else if (x <= 0.30)
		return 0.11791;
	else if (x <= 0.35)
		return 0.13683;
	else if (x <= 0.40)
		return 0.15542;
	else if (x <= 0.45)
		return 0.17365;
	else if (x <= 0.50)
		return 0.19146;
	else if (x <= 0.55)
		return 0.20884;
	else if (x <= 0.60)
		return 0.22575;
	else if (x <= 0.65)
		return 0.24215;
	else if (x <= 0.70)
		return 0.25804;
	else if (x <= 0.75)
		return 0.27337;
	else if (x <= 0.80)
		return 0.28814;
	else if (x <= 0.85)
		return 0.30234;
	else if (x <= 0.90)
		return 0.31594;
	else if (x <= 0.95)
		return 0.32894;
	else if (x <= 1.00)
		return 0.34134;
	else if (x <= 1.05)
		return 0.35314;
	else if (x <= 1.10)
		return 0.36433;
	else if (x <= 1.15)
		return 0.37493;
	else if (x <= 1.20)
		return 0.38493;
	else if (x <= 1.25)
		return 0.39435;
	else if (x <= 1.30)
		return 0.40320;
	else if (x <= 1.35)
		return 0.41149;
	else if (x <= 1.40)
		return 0.41924;
	else if (x <= 1.45)
		return 0.42647;
	else if (x <= 1.50)
		return 0.43319;
	else if (x <= 1.55)
		return 0.43943;
	else if (x <= 1.60)
		return 0.44520;
	else if (x <= 1.65)
		return 0.45053;
	else if (x <= 1.70)
		return 0.45543;
	else if (x <= 1.75)
		return 0.45994;
	else if (x <= 1.80)
		return 0.46407;
	else if (x <= 1.85)
		return 0.46784;
	else if (x <= 1.90)
		return 0.47128;
	else if (x <= 1.95)
		return 0.47441;
	else if (x <= 2.00)
		return 0.47726;
	else if (x <= 2.05)
		return 0.47982;
	else if (x <= 2.10)
		return 0.48214;
	else if (x <= 2.15)
		return 0.48422;
	else if (x <= 2.20)
		return 0.48610;
	else if (x <= 2.25)
		return 0.48778;
	else if (x <= 2.30)
		return 0.48928;
	else if (x <= 2.35)
		return 0.49061;
	else if (x <= 2.40)
		return 0.49180;
	else if (x <= 2.45)
		return 0.49286;
	else if (x <= 2.50)
		return 0.49379;
	else if (x <= 2.55)
		return 0.49461;
	else if (x <= 2.60)
		return 0.49534;
	else if (x <= 2.65)
		return 0.49597;
	else if (x <= 2.70)
		return 0.49653;
	else if (x <= 2.75)
		return 0.49702;
	else if (x <= 2.80)
		return 0.49744;
	else if (x <= 2.85)
		return 0.49781;
	else if (x <= 2.90)
		return 0.49813;
	else if (x <= 2.95)
		return 0.49841;
	else if (x <= 3.00)
		return 0.49865;
}


double UpperCase_Phi(double x)
{
	double y;
	if (x < 0)
		y = 1 - x;
	else
		y = x;
	if (y >= 4.99000000)
		return 0.9999997;
	else if (y <= -4.99000000)
		return 0.0000003;
	else
		return Normal_Distribution_CDF(y); //Phi(x);
}


double LowerCase_Phi(double x)
{
	double temp, temp1, temp2, temp3, temp4, temp5;
	temp = sqrt((2 * 3.1415));
	temp1 = 1 / temp;
	temp2 = pow(x, 2);
	temp3 = temp2 / 2;
	temp4 = -(temp3);
	temp5 = exp(temp4);
	return (temp1 * temp5);
}


double Ta(double mean1, double mean2, double teta)
{
	double temp, temp1;
	temp = (mean1 - mean2) / teta;
	//temp *= 1e9;
	if (isnan(temp))
	{
		temp1 = 0;
		return Normal_Distribution_CDF(temp1);
	}
	else if (temp < 0)
	{
		temp1 = -temp;
		if (temp1 >= 4.99000000)
			return 0.0000003;
		else
			return (1 - Normal_Distribution_CDF(temp1));
	}
	else
	{
		temp1 = temp;
		if (temp1 >= 4.99000000)
			return 0.9999997;
		else
			return Normal_Distribution_CDF(temp1);
	}
}


double Ta(double mean, double SD)
{
	double temp, temp1;
	temp = (mean + (3 * SD));
	//temp = mean;
	if (temp < 0)
		temp1 = -temp;
	else
		temp1 = temp;
	if (temp1 > 1.000000e-11)
	{
		if (temp1 > 1.000000e-10)
		{
			if (temp1 > 1.000000e-9)
			{
				if (temp1 > 1.000000e-8)
				{
					temp *= 1e8;
				}
				else
				{
					temp *= 1e9;
				}
			}
			else
			{
				temp *= 1e10;
			}
		}
		else
		{
			temp *= 1e11;
		}
	}
	else
	{
		temp *= 1e11;
	}
	if (isnan(temp))
	{
		temp1 = 0;
		return Normal_Distribution_CDF(temp1);
	}
	else if (temp < 0)
	{
		temp1 = -temp;
		if (temp1 >= 4.99000000)
			return 0.0000003;
		else
			return (1 - Normal_Distribution_CDF(temp1));
	}
	else
	{
		temp1 = temp;
		if (temp1 >= 4.99000000)
			return 0.9999997;
		else
			return Normal_Distribution_CDF(temp1);
	}
}


double Variance(vector<double> &Vth_sensitivity, double Vth_rand_sensitivity)
{
	int i;
	double temp = 0;
	for (i = 0; i < Vth_sensitivity.size(); i++)
	{
		if (Vth_sensitivity[i] != 0)
			temp += pow(Vth_sensitivity[i], 2);
	}

	return ((temp + pow(Vth_rand_sensitivity, 2)) *pow(10, -1));
}


double RO_SD1_SD2(vector<double> &a_Vth_sensitivity, vector<double> &b_Vth_sensitivity)
{
	int i;
	double temp = 0;
	for (i = 0; i < a_Vth_sensitivity.size(); i++)
		temp += (a_Vth_sensitivity[i] * b_Vth_sensitivity[i]);
	return temp;
}


double Teta(double VAR1, double VAR2, double Ro_SD1_SD2)
{
	double temp;
	temp = VAR1 + VAR2 - (2 * Ro_SD1_SD2);
	return sqrt(abs(temp));
}


#pragma once
struct Arrival
{
	vector<double> mArrival_mean;
	vector<vector<double>> mArrival_Vth_sensitivity;
	vector<double> mArrival_Vth_rand_sensitivity;
};


#pragma once
struct Sink
{
	vector<int> input_index;
	double initial_mean;
	vector<double> initial_Vth_sensitivity;
	double initial_Vth_rand_sensitivity;
	//double delay_mean;
	//vector<double> delay_Vth_sensitivity;
	//double delay_Vth_rand_sensitivity;
	double GB_mean;
	vector<double> GB_Vth_sensitivity;
	double GB_Vth_rand_sensitivity;
	double Arrival_mean;
	vector<double> Arrival_Vth_sensitivity;
	double Arrival_Vth_rand_sensitivity;
	Joint_PDF Joint;
	vector<double> Aging_Vth_sensitivity;
	vector<double> Criticality_Vth_sensitivity;
	double Aging_Vth_rand_sensitivity;
	double Criticality_Vth_rand_sensitvity;
	double RequiredArrivalTime;
	double Lifetime_RequiredArrivalTime;

	void make_input_index(vector<gate> &topological)
	{
		int i;
		for ((i = topological.size() - 1); i >= 0; i--)
		{
			if ((topological[i].output.size() == 0) && topological[i].getGateName() != "INPUT")
				input_index.push_back(i);
		}
	}


	void Guard_Band()
	{
		int i;
		GB_Vth_sensitivity.resize(16);
		GB_mean = Arrival_mean - initial_mean;
		for (i = 0; i < 16; i++)
		{
			GB_Vth_sensitivity[i] = Arrival_Vth_sensitivity[i] - initial_Vth_sensitivity[i];
		}
		GB_Vth_rand_sensitivity = Arrival_Vth_rand_sensitivity - initial_Vth_rand_sensitivity;
	}


	void Output_Delay(vector<gate> &topological, int mode = 0)
	{
		int i, j;
		double temp;
		if (mode == 0)
		{
			for (i = 0; i < input_index.size(); i++)
			{
				topological[input_index[i]].setArrival_mean(topological[input_index[i]].getDelay_mean() + topological[input_index[i]].getArrival_mean());
				for (j = 0; j < 16; j++)
					topological[input_index[i]].setArrival_Vth_sensitivity(j, topological[input_index[i]].getDelay_Vth_sensitivity(j)
					                                                                                        + topological[input_index[i]].getArrival_Vth_sensitivity(j));
				topological[input_index[i]].setArrival_Vth_rand_sensitivity(sqrt(pow(topological[input_index[i]].getDelay_Vth_rand_sensitivity(), 2)
					                                                                           + pow(topological[input_index[i]].getArrival_Vth_rand_sensitivity(), 2)));
			}
		}
		else
		{
			for (i = 0; i < input_index.size(); i++)
			{
				topological[input_index[i]].setArrival_mean(topological[input_index[i]].getDelay_mean() + topological[input_index[i]].getArrival_mean()
					                                                                                                 + topological[input_index[i]].getAging_Delay_mean());
				for (j = 0; j < 16; j++)
				{
					topological[input_index[i]].setArrival_Vth_sensitivity(j, topological[input_index[i]].getDelay_Vth_sensitivity(j)
						                     + topological[input_index[i]].getArrival_Vth_sensitivity(j) - topological[input_index[i]].getAging_Delay_Vth_sensitivity(j));
				}

				topological[input_index[i]].setArrival_Vth_rand_sensitivity(sqrt(pow(topological[input_index[i]].getDelay_Vth_rand_sensitivity()
					    - topological[input_index[i]].getAging_Delay_Vth_rand_sensitivity(), 2) + pow(topological[input_index[i]].getArrival_Vth_rand_sensitivity(), 2)));
				//topological[input_index[i]].setArrival_Vth_rand_sensitivity(sqrt(pow(topological[input_index[i]].getDelay_Vth_rand_sensitivity(), 2) 
					//+ pow(topological[input_index[i]].getArrival_Vth_rand_sensitivity(), 2) + pow(topological[input_index[i]].getAging_Delay_Vth_rand_sensitivity(), 2)));
			}
		}
	}


	void Arrival_Max(vector<gate> &topological, int mode = 0)
	{
		int i, j, k;
		double L_Phi, U_Phi, teta, ta, VAR1, VAR2, VAR_max, VAR_min, Ro_SD1_SD2, temp, temp1 = 0;

		//delay_Vth_sensitivity.resize(16);
		Arrival_Vth_sensitivity.resize(16);

		Arrival_mean = topological[input_index[0]].getArrival_mean();
		for (j = 0; j < 16; j++)
		{
			Arrival_Vth_sensitivity[j] = topological[input_index[0]].getArrival_Vth_sensitivity(j);
		}
		Arrival_Vth_rand_sensitivity = topological[input_index[0]].getArrival_Vth_rand_sensitivity();
		for (i = 1; i < input_index.size(); i++)
		{
			Arrival_mean = max(Arrival_mean, topological[input_index[i]].getArrival_mean());
			for (j = 0; j < 16; j++)
				Arrival_Vth_sensitivity[j] = max(Arrival_Vth_sensitivity[j], topological[input_index[i]].getArrival_Vth_sensitivity(j));
			Arrival_Vth_rand_sensitivity = max(Arrival_Vth_rand_sensitivity, topological[input_index[i]].getArrival_Vth_rand_sensitivity());
		}

		/*
		Arrival_mean = topological[input_index[0]].getArrival_mean();
		for (j = 0; j < 16; j++)
		{
			Arrival_Vth_sensitivity[j] = topological[input_index[0]].getArrival_Vth_sensitivity(j);
		}
		Arrival_Vth_rand_sensitivity = topological[input_index[0]].getArrival_Vth_rand_sensitivity();
		for (i = 1; i < input_index.size(); i++)
		{
			VAR1 = Variance(Arrival_Vth_sensitivity, Arrival_Vth_rand_sensitivity);
			VAR2 = Variance(topological[input_index[i]].getArrival_Vth_sensitivity(), topological[input_index[i]].getArrival_Vth_rand_sensitivity());
			Ro_SD1_SD2 = RO_SD1_SD2(Arrival_Vth_sensitivity, topological[input_index[i]].getArrival_Vth_sensitivity());
			teta = Teta(VAR1, VAR2, Ro_SD1_SD2);
			temp = (Arrival_mean - topological[input_index[i]].getArrival_mean()) / teta;
			L_Phi = LowerCase_Phi(temp);
			ta = Ta(Arrival_mean, topological[input_index[i]].getArrival_mean(), teta);
			Arrival_mean = (Arrival_mean * ta) + (topological[input_index[i]].getArrival_mean() * (1 - ta)) + (teta * L_Phi);
			VAR_max = ((VAR1 + pow(Arrival_mean, 2)) * ta) + ((VAR2 + pow(topological[input_index[i]].getArrival_mean(), 2)) * (1 - ta))
			+ ((Arrival_mean + topological[input_index[i]].getArrival_mean()) * teta * L_Phi) - pow(Arrival_mean, 2);
			for (j = 0; j < 16; j++)
			{
				Arrival_Vth_sensitivity[j] = (Arrival_Vth_sensitivity[j] * ta) + (topological[input_index[i]].getArrival_Vth_sensitivity(j) * (1 - ta));
				temp1 += pow(Arrival_Vth_sensitivity[j], 2);
			}
			Arrival_Vth_rand_sensitivity = sqrt(VAR_max - temp1);
			temp1 = 0;
		}
		*/
	}


	void Joint_PDF_Computation(vector<gate> &topological)
	{
		int j;
		
		Aging_Vth_sensitivity.resize(16);
		Criticality_Vth_sensitivity.resize(16);
		for (j = 0; j < 16; j++)
		{
			Aging_Vth_sensitivity[j] = Arrival_Vth_sensitivity[j];// -initial_Vth_sensitivity[j];
			Criticality_Vth_sensitivity[j] = Arrival_Vth_sensitivity[j] - topological[input_index[0]].getRAT_Vth_sensitivity(j);
		}
		Aging_Vth_rand_sensitivity = sqrt(pow(Arrival_Vth_rand_sensitivity, 2));// +pow(initial_Vth_rand_sensitivity, 2));
		Criticality_Vth_rand_sensitvity = sqrt(pow(Arrival_Vth_rand_sensitivity, 2) + pow(topological[input_index[0]].getRAT_Vth_rand_sensitivity(), 2));
		Joint.Aging_mean = Arrival_mean - Lifetime_RequiredArrivalTime;// initial_mean;
		Joint.Criticality_mean = Arrival_mean - topological[input_index[0]].getRAT_mean();
		Joint.Aging_Variance = Variance(Aging_Vth_sensitivity, Aging_Vth_rand_sensitivity);
		Joint.Criticality_Variance = Variance(Criticality_Vth_sensitivity, Criticality_Vth_rand_sensitvity);
		Joint.Correlation_coefficient = RO_SD1_SD2(Aging_Vth_sensitivity, Criticality_Vth_sensitivity);
	}
};


void SSTA_ArrivalTime(vector<gate> &topological, Sink *sink, int mode = 0)
{
	int i, j, k;
	double L_Phi, teta, ta, VAR1, VAR2, VAR_max, Ro_SD1_SD2;
	double temp, temp1 = 0, temp2;
	vector<Arrival> AT;
	AT.resize(topological.size());

	for (i = 0; i < topological.size(); i++)
	{
		if (topological[i].getGateName() != "INPUT")
		{
			// max

			topological[i].setArrival_mean(AT[i].mArrival_mean[0]);
			for (k = 0; k < 16; k++)
			{
				topological[i].setArrival_Vth_sensitivity(k, AT[i].mArrival_Vth_sensitivity[0][k]);
			}
			topological[i].setArrival_Vth_rand_sensitivity(AT[i].mArrival_Vth_rand_sensitivity[0]);
			for (j = 1; j < AT[i].mArrival_mean.size(); j++)
			{
				topological[i].setArrival_mean(max(topological[i].getArrival_mean(), AT[i].mArrival_mean[j]));
				for (k = 0; k < 16; k++)
					topological[i].setArrival_Vth_sensitivity(k, max(topological[i].getArrival_Vth_sensitivity(k), AT[i].mArrival_Vth_sensitivity[j][k]));
				topological[i].setArrival_Vth_rand_sensitivity(max(topological[i].getArrival_Vth_rand_sensitivity(), AT[i].mArrival_Vth_rand_sensitivity[j]));
			}

			/*
			topological[i].setArrival_mean(AT[i].mArrival_mean[0]);
			for (k = 0; k < 16; k++)
			{
				topological[i].setArrival_Vth_sensitivity(k, AT[i].mArrival_Vth_sensitivity[0][k]);
			}
			topological[i].setArrival_Vth_rand_sensitivity(AT[i].mArrival_Vth_rand_sensitivity[0]);
			for (j = 1; j < AT[i].mArrival_mean.size(); j++)
			{
				VAR1 = Variance(topological[i].getArrival_Vth_sensitivity(), topological[i].getArrival_Vth_rand_sensitivity());
				VAR2 = Variance(AT[i].mArrival_Vth_sensitivity[j], AT[i].mArrival_Vth_rand_sensitivity[j]);
				Ro_SD1_SD2 = RO_SD1_SD2(topological[i].getArrival_Vth_sensitivity(), AT[i].mArrival_Vth_sensitivity[j]);
				teta = Teta(VAR1, VAR2, Ro_SD1_SD2);
				temp = (topological[i].getArrival_mean() - AT[i].mArrival_mean[j]) / teta;
				if (isnan(temp))
					temp = 0;
				L_Phi = LowerCase_Phi(temp);
				ta = Ta(topological[i].getArrival_mean(), AT[i].mArrival_mean[j], teta);
				topological[i].setArrival_mean((topological[i].getArrival_mean() * ta) + (AT[i].mArrival_mean[j] * (1 - ta)) + (teta * L_Phi));
				VAR_max = ((VAR1 + pow(topological[i].getArrival_mean(), 2)) * ta) + ((VAR2 + pow(AT[i].mArrival_mean[j], 2)) * (1 - ta))
					+ ((topological[i].getArrival_mean() + AT[i].mArrival_mean[j]) * teta * L_Phi) - pow(topological[i].getArrival_mean(), 2);
				for (k = 0; k < 16; k++)
				{
					topological[i].setArrival_Vth_sensitivity(k, (topological[i].getArrival_Vth_sensitivity(k) * ta) + (AT[i].mArrival_Vth_sensitivity[j][k] * (1 - ta)));
					temp1 += pow(topological[i].getArrival_Vth_sensitivity(k), 2);
				}
				topological[i].setArrival_Vth_rand_sensitivity(sqrt(VAR_max - temp1));
				temp1 = 0;
			}
			*/
			// sum
			if (mode == 0)
			{
				for (j = 0; j < topological[i].output.size(); j++)
				{
					AT[topological[i].output[j]].mArrival_mean.push_back(topological[i].getDelay_mean() + topological[i].getArrival_mean());
					AT[topological[i].output[j]].mArrival_Vth_sensitivity.resize(AT[topological[i].output[j]].mArrival_Vth_sensitivity.size() + 1);
					for (k = 0; k < 16; k++)
					{
						AT[topological[i].output[j]].mArrival_Vth_sensitivity[AT[topological[i].output[j]].mArrival_Vth_sensitivity.size() - 1].push_back(
							                                                   topological[i].getDelay_Vth_sensitivity(k) + topological[i].getArrival_Vth_sensitivity(k));
					}
					AT[topological[i].output[j]].mArrival_Vth_rand_sensitivity.push_back(sqrt(pow(topological[i].getDelay_Vth_rand_sensitivity(), 2)
						                                                                                     + pow(topological[i].getArrival_Vth_rand_sensitivity(), 2)));
				}
			}
			else
			{
				for (j = 0; j < topological[i].output.size(); j++)
				{
					AT[topological[i].output[j]].mArrival_mean.push_back(topological[i].getDelay_mean() + topological[i].getArrival_mean() 
						                                                                                                          + topological[i].getAging_Delay_mean());
					AT[topological[i].output[j]].mArrival_Vth_sensitivity.resize(AT[topological[i].output[j]].mArrival_Vth_sensitivity.size() + 1);
					for (k = 0; k < 16; k++)
					{
						AT[topological[i].output[j]].mArrival_Vth_sensitivity[AT[topological[i].output[j]].mArrival_Vth_sensitivity.size() - 1].push_back(
							topological[i].getDelay_Vth_sensitivity(k) + topological[i].getArrival_Vth_sensitivity(k) - topological[i].getAging_Delay_Vth_sensitivity(k));
					}

					AT[topological[i].output[j]].mArrival_Vth_rand_sensitivity.push_back(sqrt(pow(topological[i].getDelay_Vth_rand_sensitivity()
						                          - topological[i].getAging_Delay_Vth_rand_sensitivity(), 2) + pow(topological[i].getArrival_Vth_rand_sensitivity(), 2)));
					//AT[topological[i].output[j]].mArrival_Vth_rand_sensitivity.push_back(sqrt(pow(topological[i].getDelay_Vth_rand_sensitivity(), 2) 
						                      //+ pow(topological[i].getArrival_Vth_rand_sensitivity(), 2) + pow(topological[i].getAging_Delay_Vth_rand_sensitivity(), 2))); 
				}
			}
		}
		else
		{
			for (j = 0; j < topological[i].output.size(); j++)
			{
				AT[topological[i].output[j]].mArrival_mean.push_back(0);
				AT[topological[i].output[j]].mArrival_Vth_sensitivity.resize(AT[topological[i].output[j]].mArrival_Vth_sensitivity.size() + 1);
				for (k = 0; k < 16; k++)
				{
					AT[topological[i].output[j]].mArrival_Vth_sensitivity[AT[topological[i].output[j]].mArrival_Vth_sensitivity.size() - 1].push_back(0);
				}
				AT[topological[i].output[j]].mArrival_Vth_rand_sensitivity.push_back(0);
			}
		}
	}
}


void SSTA_ReqiredArrivalTime(vector<gate> &topological, Sink *sink, int mode = 0)
{
	int i, j, k;
	double L_Phi, teta, ta, VAR1, VAR2, VAR_min, Ro_SD1_SD2;
	double temp, temp1 = 0;
	vector<Arrival> AT;
	AT.resize(topological.size());

	for ((i = topological.size() - 1); i >= 0; i--)
	{
		if ((topological[i].output.size() != 0) && (topological[i].getGateName() != "INPUT"))
		{
			// min

			topological[i].setRAT_mean(AT[i].mArrival_mean[0]);
			for (k = 0; k < 16; k++)
			{
				topological[i].setRAT_Vth_sensitivity(k, AT[i].mArrival_Vth_sensitivity[0][k]);
			}
			topological[i].setRAT_Vth_rand_sensitivity(AT[i].mArrival_Vth_rand_sensitivity[0]);
			for (j = 1; j < AT[i].mArrival_mean.size(); j++)
			{
				topological[i].setRAT_mean(min(topological[i].getRAT_mean(), AT[i].mArrival_mean[j]));
				for (k = 0; k < 16; k++)
					topological[i].setRAT_Vth_sensitivity(k, min(topological[i].getRAT_Vth_sensitivity(k), AT[i].mArrival_Vth_sensitivity[j][k]));
				topological[i].setRAT_Vth_rand_sensitivity(min(topological[i].getRAT_Vth_rand_sensitivity(), AT[i].mArrival_Vth_rand_sensitivity[j]));
			}

			/*
			topological[i].setRAT_mean(AT[i].mArrival_mean[0]);
			for (k = 0; k < 16; k++)
			{
				topological[i].setRAT_Vth_sensitivity(k, AT[i].mArrival_Vth_sensitivity[0][k]);
			}
			topological[i].setRAT_Vth_rand_sensitivity(AT[i].mArrival_Vth_rand_sensitivity[0]);
			for (j = 1; j < AT[i].mArrival_mean.size(); j++)
			{
				VAR1 = Variance(topological[i].getRAT_Vth_sensitivity(), topological[i].getRAT_Vth_rand_sensitivity());
				VAR2 = Variance(AT[i].mArrival_Vth_sensitivity[j], AT[i].mArrival_Vth_rand_sensitivity[j]);
				Ro_SD1_SD2 = RO_SD1_SD2(topological[i].getRAT_Vth_sensitivity(), AT[i].mArrival_Vth_sensitivity[j]);
				teta = Teta(VAR1, VAR2, Ro_SD1_SD2);
				temp = (AT[i].mArrival_mean[j] - topological[i].getRAT_mean()) / teta;
				if (isnan(temp))
					temp = 0;
				L_Phi = LowerCase_Phi(temp);
				ta = Ta(AT[i].mArrival_mean[j], topological[i].getRAT_mean(), teta);
				topological[i].setRAT_mean((topological[i].getRAT_mean() * ta) + (AT[i].mArrival_mean[j] * (1 - ta)) - (teta * L_Phi));
				VAR_min = ((VAR1 + pow(topological[i].getRAT_mean(), 2)) * ta) + ((VAR2 + pow(AT[i].mArrival_mean[j], 2)) * (1 - ta))
					- ((topological[i].getRAT_mean() + AT[i].mArrival_mean[j]) * teta * L_Phi) - pow(topological[i].getRAT_mean(), 2);
				for (k = 0; k < 16; k++)
				{
					topological[i].setRAT_Vth_sensitivity(k, abs((topological[i].getRAT_Vth_sensitivity(k) * ta) + (AT[i].mArrival_Vth_sensitivity[j][k] * (1 - ta))));
					temp1 += pow(topological[i].getRAT_Vth_sensitivity(k), 2);
				}

				topological[i].setRAT_Vth_rand_sensitivity(sqrt(abs((VAR_min - temp1))));
				temp1 = 0;
			}
			*/

			//sub         //*************************************** correct
			if (mode == 0)
			{
				for (j = 0; j < topological[i].getIndex_Inputs().size(); j++)
				{
					AT[topological[i].getIndex_Inputs()[j]].mArrival_mean.push_back(topological[i].getRAT_mean() - topological[i].getDelay_mean());
					AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_sensitivity.resize(AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_sensitivity.size() + 1);
					for (k = 0; k < 16; k++)
					{
						AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_sensitivity[AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_sensitivity.size() - 1].push_back(
							topological[i].getRAT_Vth_sensitivity(k) - topological[i].getDelay_Vth_sensitivity(k));
					}
					AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_rand_sensitivity.push_back(
						                                                         sqrt(pow(topological[i].getDelay_mean(), 2) + pow(topological[i].getArrival_mean(), 2)));
				}
			}
			else
			{
				for (j = 0; j < topological[i].getIndex_Inputs().size(); j++)
				{
					AT[topological[i].getIndex_Inputs()[j]].mArrival_mean.push_back(
						                                          topological[i].getRAT_mean() - (topological[i].getDelay_mean() + topological[i].getAging_Delay_mean()));
					AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_sensitivity.resize(AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_sensitivity.size() + 1);
					for (k = 0; k < 16; k++)
					{
						AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_sensitivity[AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_sensitivity.size() - 1].push_back(
							topological[i].getRAT_Vth_sensitivity(k) - (topological[i].getDelay_Vth_sensitivity(k) - topological[i].getAging_Delay_Vth_sensitivity(k)));
					}

					AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_rand_sensitivity.push_back(
						sqrt(pow(topological[i].getDelay_mean() - topological[i].getAging_Delay_mean(), 2) + pow(topological[i].getRAT_mean(), 2)));
					//AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_rand_sensitivity.push_back(
						//sqrt(pow(topological[i].getDelay_mean(), 2) + pow(topological[i].getRAT_mean(), 2) + pow(topological[i].getAging_Delay_mean(), 2)));
				}
			}
		}
		else if ((topological[i].getGateName() == "INPUT") || ((topological[i].output.size() == 0) && topological[i].getGateName() == "INPUT"))
			continue;
		else
		{
			for (j = 0; j < topological[i].getIndex_Inputs().size(); j++)
			{
				AT[topological[i].getIndex_Inputs()[j]].mArrival_mean.push_back(topological[i].getRAT_mean());
				AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_sensitivity.resize(AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_sensitivity.size() + 1);
				for (k = 0; k < 16; k++)
				{
					AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_sensitivity[AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_sensitivity.size() - 1]
						                                                                                             .push_back(topological[i].getRAT_Vth_sensitivity(k));
				}
				AT[topological[i].getIndex_Inputs()[j]].mArrival_Vth_rand_sensitivity.push_back(topological[i].getRAT_Vth_rand_sensitivity());
			}
		}
	}
}


void SSTA(vector<gate> &topological, Sink *sink, double initial_constr, double lifetime_constr, mode Mode = initial)
{
	int i, j;
	double VAR1, VAR2, Ro_SD1_SD2, teta;
	if (Mode == initial)
	{
		SSTA_ArrivalTime(topological, sink, 0);// Mode);
		sink->make_input_index(topological);
		//**** calculate RAT
		sink->Output_Delay(topological, 0);// Mode);
		sink->Arrival_Max(topological, 0);// Mode);
		//sink->Initial_Constraint(initial_constr * (sink->Arrival_mean + (3 * sqrt(Variance(sink->Arrival_Vth_sensitivity, sink->Arrival_Vth_rand_sensitivity)))));
		//sink->Lifetime_Constraint(constr);
		sink->RequiredArrivalTime = initial_constr * (sink->Arrival_mean + (3 * sqrt(Variance(sink->Arrival_Vth_sensitivity, sink->Arrival_Vth_rand_sensitivity))));
		sink->Lifetime_RequiredArrivalTime = lifetime_constr * sink->RequiredArrivalTime;
		/*VAR1 = Variance(sink->initial_delay_Vth_sensitivity, sink->initial_delay_Vth_rand_sensitivity);
		VAR2 = Variance(sink->Arrival_Vth_sensitivity, sink->Arrival_Vth_rand_sensitivity);
		Ro_SD1_SD2 = RO_SD1_SD2(sink->initial_delay_Vth_sensitivity, sink->Arrival_Vth_sensitivity);
		teta = Teta(VAR1, VAR2, Ro_SD1_SD2);
		*Constraint = Ta(sink->initial_delay_mean, sink->Arrival_mean, teta);*/
		sink->initial_mean = sink->Arrival_mean;
		for (i = 0; i < 16; i++)
			sink->initial_Vth_sensitivity.push_back(sink->Arrival_Vth_sensitivity[i]);
		sink->initial_Vth_rand_sensitivity = sink->Arrival_Vth_rand_sensitivity;
		for (i = 0; i < sink->input_index.size(); i++)
		{
			topological[sink->input_index[i]].setRAT_mean(initial_constr * (sink->Arrival_mean + (3 * sqrt(Variance(sink->Arrival_Vth_sensitivity, sink->Arrival_Vth_rand_sensitivity)))));
			
			for (j = 0; j < 16; j++)
			{
				topological[sink->input_index[i]].setRAT_Vth_sensitivity(j, 0);
			}
			topological[sink->input_index[i]].setRAT_Vth_rand_sensitivity(0);
			
		}
		SSTA_ReqiredArrivalTime(topological, sink, 0);// Mode);
	}
	else if (Mode == initial_optimization)
	{
		SSTA_ArrivalTime(topological, sink, 0);// Mode);
		sink->Output_Delay(topological, 0);// Mode);
		sink->Arrival_Max(topological, 0);// Mode);
		SSTA_ReqiredArrivalTime(topological, sink, 0);// Mode);
		sink->initial_mean = sink->Arrival_mean;
		for (i = 0; i < 16; i++)
			sink->initial_Vth_sensitivity.push_back(sink->Arrival_Vth_sensitivity[i]);
		sink->initial_Vth_rand_sensitivity = sink->Arrival_Vth_rand_sensitivity;
	}
	else if (Mode == lifetime)
	{
		SSTA_ArrivalTime(topological, sink, 1);// Mode);
		sink->Output_Delay(topological, 1);// Mode);
		sink->Arrival_Max(topological, 1);// Mode);
		for (i = 0; i < sink->input_index.size(); i++)
		{
			topological[sink->input_index[i]].setRAT_mean(sink->Lifetime_RequiredArrivalTime);

			for (j = 0; j < 16; j++)
			{
				topological[sink->input_index[i]].setRAT_Vth_sensitivity(j, 0);
			}
			topological[sink->input_index[i]].setRAT_Vth_rand_sensitivity(0);

		}
		sink->Guard_Band();
		SSTA_ReqiredArrivalTime(topological, sink, 1);// Mode);
	}
	else
	{
		SSTA_ArrivalTime(topological, sink, 1);// Mode);
		sink->Output_Delay(topological, 1);// Mode);
		sink->Arrival_Max(topological, 1);// Mode);
		SSTA_ReqiredArrivalTime(topological, sink, 1);// Mode);
	}
}


void Criticality_Computation(vector<gate> &topological, Sink *sink)
{
	int i, j;
	double VAR1, VAR2, Ro_SD1_SD2, teta, U_Phi;
	double mean, Vth_rand_sensitivity, dVth_sensitivity;
	vector<double> Vth_sensitivity;
	for (i = 0; i < topological.size(); i++)
	{
		if (topological[i].getGateName() != "INPUT")
		{
			Vth_sensitivity.clear();
			Vth_sensitivity.resize(16);
			mean = topological[i].getArrival_mean() - topological[i].getRAT_mean();
			for (j = 0; j < 16; j++)
			{
				Vth_sensitivity[j] = topological[i].getArrival_Vth_sensitivity(j) - topological[i].getRAT_Vth_sensitivity(j);
			}
			//****************************************************
			Vth_rand_sensitivity = sqrtl(pow(topological[i].getArrival_Vth_rand_sensitivity(), 2) + pow(topological[i].getRAT_Vth_rand_sensitivity(), 2));
			//****************************************************
			//Vth_rand_sensitivity = topological[i].getArrival_Vth_rand_sensitivity() - topological[i].getRAT_Vth_rand_sensitivity();
			VAR1 = Variance(Vth_sensitivity, Vth_rand_sensitivity);
			VAR2 = Variance(sink->Arrival_Vth_sensitivity, sink->Arrival_Vth_rand_sensitivity);
			Ro_SD1_SD2 = RO_SD1_SD2(Vth_sensitivity, sink->Arrival_Vth_sensitivity);
			teta = Teta(VAR1, VAR2, Ro_SD1_SD2);
			topological[i].setCriticality(Ta(mean, sink->Arrival_mean, teta));
		}
	}
}


void Average_Criticality_Computation(vector<gate> &topological)
{
	int i, j;
	double temp;

	for (i = topological.size() - 1; i >= 0; i--)
	{
		if (topological[i].output.size() == 0)
			topological[i].setAverage_Criticality(topological[i].getCriticality());
	}

	for (i = topological.size() - 1; i >= 0; i--)
	{
		temp = 0;
		if ((topological[i].output.size() != 0) && (topological[i].getGateName() != "INPUT"))
		{
			for (j = 0; j < topological[i].output.size(); j++)
				temp += topological[topological[i].output[j]].getAverage_Criticality();
			topological[i].setAverage_Criticality(temp / topological[i].output.size());
		}
	}
}


void Rank_Criticality(vector<gate> &topological, vector<Critical> &Criticality)
{
	int i;
	Criticality.resize(topological.size());
	for (i = 0; i < topological.size(); i++)
	{
		Criticality[i].criticality = topological[i].getCriticality();
		Criticality[i].index = i;
	}

	stable_sort(Criticality.begin(), Criticality.end(), [](const Critical &a, const Critical &b)  {	return (a.criticality < b.criticality);} );
	reverse(Criticality.begin(), Criticality.end());
}


void Rank_Average_Criticality(vector<gate> &topological, vector<Critical> &Criticality)
{
	int i;
	Criticality.resize(topological.size());
	for (i = 0; i < topological.size(); i++)
	{
		Criticality[i].criticality = topological[i].getAverage_Criticality();
		Criticality[i].index = i;
	}

	stable_sort(Criticality.begin(), Criticality.end(), [](const Critical &a, const Critical &b)  {	return (a.criticality < b.criticality); });
	reverse(Criticality.begin(), Criticality.end());
}


void Aging_Computation(vector<gate> &topological, Sink *sink)
{
	int i;
	double VAR_aging, VAR_GB;
	double Ro_SD1_SD2, teta;
	//double temp;

	VAR_GB = Variance(sink->GB_Vth_sensitivity, sink->GB_Vth_rand_sensitivity);
	for (i = 0; i < topological.size(); i++)
	{
		if (topological[i].getGateName() != "INPUT")
		{
			VAR_aging = Variance(topological[i].getAging_Delay_Vth_sensitivity(), topological[i].getAging_Delay_Vth_rand_sensitivity());
			Ro_SD1_SD2 = RO_SD1_SD2(topological[i].getAging_Delay_Vth_sensitivity(), sink->GB_Vth_sensitivity);
			teta = Teta(VAR_aging, VAR_GB, Ro_SD1_SD2);
			//temp = (topological[i].getAging_Delay_mean() - sink->GB_mean) / teta;
			topological[i].setAging(Ta(topological[i].getAging_Delay_mean(), sink->GB_mean, teta));
		}
	}
}


void Pareto_Ranking(vector<gate> &topological, Sink *sink, vector<Age_Critical> &Aging_Critical)
{
	int i, j, k;
	for (i = 0; i < topological.size(); i++)
	{
		k = 0;
		if (topological[i].getGateName() != "INPUT")
		{
			for (j = 0; j < topological.size(); j++)
			{
				if ((topological[j].getGateName() != "INPUT") && (i != j))
				{
					if (topological[i].getCriticality() > topological[j].getCriticality() || (topological[i].getAging() > topological[j].getAging()))
						k++;
				}
			}
			topological[i].setRank(k);
			Aging_Critical.resize(Aging_Critical.size() + 1);
			Aging_Critical[Aging_Critical.size() - 1].index = i;
			Aging_Critical[Aging_Critical.size() - 1].Pareto_index = k;
		}
	}
	stable_sort(Aging_Critical.begin(), Aging_Critical.end(), [](const Age_Critical &a, const Age_Critical &b)  {	return (a.Pareto_index < b.Pareto_index); });
	reverse(Aging_Critical.begin(), Aging_Critical.end());
}


void Pareto_Ranking_Multiply(vector<gate> &topological, Sink *sink, vector<Critical> &Criticality)
{
	int i;
	Criticality.resize(topological.size());
	for (i = 0; i < topological.size(); i++)
	{
		if (topological[i].getGateName() != "INPUT")
		{
			Criticality[i].criticality = topological[i].getCriticality() * topological[i].getAging();
			Criticality[i].index = i;
		}
		else
		{
			Criticality[i].criticality = 0;
			Criticality[i].index = i;
		}
	}
	stable_sort(Criticality.begin(), Criticality.end(), [](const Critical &a, const Critical &b)  {	return (a.criticality < b.criticality); });
	reverse(Criticality.begin(), Criticality.end());
}


void Criticality_PDF(vector<gate> &topological, int i)
{
	int j;
	double VAR1, VAR2, Ro_SD1_SD2, teta, U_Phi;
	double mean, Vth_rand_sensitivity, dVth_sensitivity;
	vector<double> Vth_sensitivity;

	topological[i].setCriticality_mean(topological[i].getArrival_mean() - topological[i].getRAT_mean());
	for (j = 0; j < 16; j++)
	{
		topological[i].setCriticality_Vth_sensitivity(j, topological[i].getArrival_Vth_sensitivity(j) - topological[i].getRAT_Vth_sensitivity(j));
	}
	topological[i].setCriticality_Vth_rand_sensitivity(sqrt(pow(topological[i].getArrival_Vth_rand_sensitivity(), 2)
		 	                                                                                                     + pow(topological[i].getRAT_Vth_rand_sensitivity(), 2)));
		
	/*
	int i, j;
	double VAR1, VAR2, Ro_SD1_SD2, teta, U_Phi;
	double mean, Vth_rand_sensitivity, dVth_sensitivity;
	vector<double> Vth_sensitivity;
	for (i = 0; i < topological.size(); i++)
	{
		if (topological[i].getGateName() != "INPUT")
		{
			
			topological[i].setCriticality_mean = topological[i].getArrival_mean() - topological[i].getRAT_mean();
			for (j = 0; j < 16; j++)
			{
				topological[i].setCriticality_Vth_sensitivity(j, topological[i].getArrival_Vth_sensitivity(j) - topological[i].getRAT_Vth_sensitivity(j));
			}
			topological[i].setCriticality_Vth_rand_sensitivity(sqrt(pow(topological[i].getArrival_Vth_rand_sensitivity(), 2) 
				                                                                                                 + pow(topological[i].getRAT_Vth_rand_sensitivity(), 2)));
		}
	}
	*/
}


void Join_PDF_Computation(vector<gate> &topological, int i)
{
	double VAR1, VAR2, Ro_SD1_SD2, teta, U_Phi;
	double mean, Vth_rand_sensitivity, dVth_sensitivity;
	Joint_PDF Joint;
	vector<double> Vth_sensitivity;
	
	VAR1 = Variance(topological[i].getAging_Delay_Vth_sensitivity(), topological[i].getAging_Delay_Vth_rand_sensitivity());
	VAR2 = Variance(topological[i].getCriticality_Vth_sensitivity(), topological[i].getCriticality_Vth_rand_sensitivity());
	Ro_SD1_SD2 = RO_SD1_SD2(topological[i].getAging_Delay_Vth_sensitivity(), topological[i].getCriticality_Vth_sensitivity());
	Joint.Aging_mean = topological[i].getAging_Delay_mean();
	Joint.Criticality_mean = topological[i].getCriticality_mean();
	Joint.Aging_Variance = VAR1;
	Joint.Criticality_Variance = VAR2;
	Joint.Correlation_coefficient = Ro_SD1_SD2;// (Ro_SD1_SD2 / (sqrt(VAR1) * sqrt(VAR2)));
	topological[i].setJoint(Joint);
	/*
	int i, j;
	double VAR1, VAR2, Ro_SD1_SD2, teta, U_Phi;
	double mean, Vth_rand_sensitivity, dVth_sensitivity;
	Joint_PDF Joint;
	vector<double> Vth_sensitivity;
	for (i = 0; i < topological.size(); i++)
	{
		if (topological[i].getGateName() != "INPUT")
		{
			VAR1 = Variance(topological[i].getAging_Delay_Vth_sensitivity(), topological[i].getAging_Delay_Vth_rand_sensitivity());
			VAR2 = Variance(topological[i].getCriticality_Vth_sensitivity(), topological[i].getCriticality_Vth_rand_sensitivity());
			Ro_SD1_SD2 = RO_SD1_SD2(topological[i].getAging_Delay_Vth_sensitivity(), topological[i].getCriticality_Vth_sensitivity());
			Joint.Aging_mean = topological[i].getAging_Delay_mean();
			Joint.Criticality_mean = topological[i].getCriticality_mean();
			Joint.Aging_Variance = VAR1;
			Joint.Criticality_Variance = VAR2;
			Joint.Correlation_coefficient = (Ro_SD1_SD2 / (sqrt(VAR1) * sqrt(VAR2)));
			topological[i].setJoint(Joint);
		}
	}
	*/
}


Joint_PDF Joint_PDF_SUB(Joint_PDF i, Joint_PDF s, vector<double> i_Aging_Vth_sensitivity, vector<double> i_Criticality_Vth_sensitvity, vector<double> s_Aging_Vth_sensitivity,
	                                                                                                                vector<double> s_Criticality_Vth_sensitivity)
{
	int j;
	double cov1, cov2;
	Joint_PDF Joint;
	vector<double> Aging_Vth_sensitivity, Criticality_Vth_sensitivity;
	cov1 = RO_SD1_SD2(i_Aging_Vth_sensitivity, s_Aging_Vth_sensitivity);
	cov2 = RO_SD1_SD2(i_Criticality_Vth_sensitvity, s_Criticality_Vth_sensitivity);
	Joint.Aging_mean = i.Aging_mean - s.Aging_mean;
	Joint.Criticality_mean = i.Criticality_mean - s.Criticality_mean;
	Joint.Aging_Variance = i.Aging_Variance + s.Aging_Variance - (2 * cov1);
	Joint.Criticality_Variance = i.Criticality_Variance + s.Criticality_Variance - (2 * cov2);
	Aging_Vth_sensitivity.resize(16);
	Criticality_Vth_sensitivity.resize(16);
	for (j = 0; j < 16; j++)
	{
		Aging_Vth_sensitivity[j] = i_Aging_Vth_sensitivity[j] - s_Aging_Vth_sensitivity[j];
		Criticality_Vth_sensitivity[j] = i_Criticality_Vth_sensitvity[j] - s_Criticality_Vth_sensitivity[j];
	}
	Joint.Correlation_coefficient = RO_SD1_SD2(Aging_Vth_sensitivity, Criticality_Vth_sensitivity);
	return Joint;
}


void Critical_Aging_Computation(vector<gate> &topological, Sink *sink)
{
	int i;
	Engine *ep;
	Joint_PDF Joint;
	for (i = 0; i < topological.size(); i++)
	{
		if (topological[i].getGateName() != "INPUT")
		{
			Criticality_PDF(topological, i);
			Join_PDF_Computation(topological, i);
		}
	}
	sink->Joint_PDF_Computation(topological);

	if (!(ep = engOpen(NULL)))
	{
		cout << "error" << endl;
	}

	for (i = 0; i < topological.size(); i++)
	{
		if (topological[i].getGateName() != "INPUT")
		{
			Joint = Joint_PDF_SUB(topological[i].getJoint_PDF(), sink->Joint, topological[i].getAging_Delay_Vth_sensitivity(), topological[i].getCriticality_Vth_sensitivity(),
				sink->Aging_Vth_sensitivity, sink->Criticality_Vth_sensitivity);
			double xArray[2] = { 0, 0 };
			xArray[0] = 0;
			xArray[1] = 0;
			mxArray* x = mxCreateDoubleMatrix(2, 1, mxREAL);
			memcpy((void *)mxGetPr(x), (void *)xArray, sizeof(xArray));
			engPutVariable(ep, "x", x);

			double muArray[2];
			muArray[0] = Joint.Aging_mean;
			muArray[1] = Joint.Criticality_mean;
			mxArray* mu = mxCreateDoubleMatrix(2, 1, mxREAL);
			memcpy((void *)mxGetPr(mu), (void *)muArray, sizeof(muArray));
			engPutVariable(ep, "mu", mu);

			double sigmaArray[2][2];
			sigmaArray[0][0] = Joint.Aging_Variance;
			sigmaArray[0][1] = Joint.Correlation_coefficient;
			sigmaArray[1][0] = Joint.Correlation_coefficient;
			sigmaArray[1][1] = Joint.Criticality_Variance;
			mxArray* sigma = mxCreateDoubleMatrix(2, 2, mxREAL);
			memcpy((void *)mxGetPr(sigma), (void *)sigmaArray, sizeof(sigmaArray));
			engPutVariable(ep, "sigma", sigma);

			engEvalString(ep, "p=mvncdf(x,mu,sigma);");

			double *cresult;
			mxArray *mresult;
			mresult = engGetVariable(ep, "p");
			cresult = mxGetPr(mresult);
			topological[i].setCritical_Aging(1 - (*cresult));
		}
		else
			topological[i].setCritical_Aging(0);
	}
	//engClose(ep);
	/*
	double xArray[2] = { 0, 0 };
	mxArray* x = mxCreateDoubleMatrix(2, 1, mxREAL);
	memcpy((void *)mxGetPr(x), (void *)xArray, sizeof(xArray));
	engPutVariable(ep, "x", x);

	double muArray[2];
	mxArray* mu = mxCreateDoubleMatrix(2, 1, mxREAL);
	memcpy((void *)mxGetPr(mu), (void *)muArray, sizeof(muArray));
	engPutVariable(ep, "mu", mu);
	
	double sigmaArray[2][2] = { 1, 0, 0, 1 };
	mxArray* sigma = mxCreateDoubleMatrix(2, 2, mxREAL);
	memcpy((void *)mxGetPr(sigma), (void *)sigmaArray, sizeof(sigmaArray));
	engPutVariable(ep, "sigma", sigma);
	

	//engEvalString(ep, "x=[0 0];");
	//engEvalString(ep, "mu=[0 0];");
	//engEvalString(ep, "sigma=[1 0; 0 1];");
	engEvalString(ep, "p=mvncdf(x,mu,sigma);");
	engEvalString(ep, "plot(x,mu);");
	double *cresult;
	mxArray *mresult;
	mresult = engGetVariable(ep, "p");
	cresult = mxGetPr(mresult);
	*/
}


void Rank_Critical_Aging(vector<gate> &topological, vector<Critical> &Criticality)
{
	int i;
	Criticality.resize(topological.size());
	for (i = 0; i < topological.size(); i++)
	{
		Criticality[i].criticality = topological[i].getCritical_Aging();
		Criticality[i].index = i;
	}

	stable_sort(Criticality.begin(), Criticality.end(), [](const Critical &a, const Critical &b)  {	return (a.criticality < b.criticality); });
	reverse(Criticality.begin(), Criticality.end());
}


void Spatial_Correlation(vector<gate> &topological, vector<int> input)
{
	int i, j = 1, k = 5, z = 0, m = 0;
	double temp = ceil((topological.size() - input.size()) / 4.0);
	double temp1 = ceil((topological.size() - input.size()) / 16.0);
	for (i = 0; i < topological.size(); i++)
	{
		if (topological[i].getGateName() != "INPUT")
		{
			topological[i].setDelay_mean();
			topological[i].setAging_Delay_mean(10 * 365 * 24 * 60 * 60);
			topological[i].setDelay_Vth_rand_sensitivity();
			topological[i].setAging_Delay_Vth_rand_sensitivity(10 * 365 * 24 * 60 * 60);
			topological[i].setDelay_Vth_sensitivity(0);
			topological[i].setAging_Delay_Vth_sensitivity(0, 10 * 365 * 24 * 60 * 60);
			topological[i].setDelay_Vth_sensitivity(j);
			topological[i].setAging_Delay_Vth_sensitivity(j, 10 * 365 * 24 * 60 * 60);
			z++;
			if (z == temp)
			{
				z = 0;
				j++;
			}
			topological[i].setDelay_Vth_sensitivity(k);
			topological[i].setAging_Delay_Vth_sensitivity(k, 10 * 365 * 24 * 60 * 60);
			m++;
			if (m == temp1)
			{
				m = 0;
				if (k != 15)
					k++;
			}
		}
	}
}


double Yield_Computation(double Delay_mean, vector<double> &Delay_Vth, double Delay_Vth_rand,double Constraint_mean, vector<double> &Constraint_Vth,
	                                                                                                                                          double Constraiant_Vth_rand)
{
	int i;
	double Vth_rand_sensitivity, VAR, VAR1, VAR2, Ro_SD1_SD2;
	vector<double> Vth_sensitivity;

	Vth_sensitivity.resize(16);

	for (i = 0; i < 16; i++)
		Vth_sensitivity[i] = Constraint_Vth[i] - Delay_Vth[i];
	//Vth_rand_sensitivity = Constraiant_Vth_rand - Delay_Vth_rand;
	VAR1 = Variance(Delay_Vth, Delay_Vth_rand);
	VAR2 = Variance(Constraint_Vth, Constraiant_Vth_rand);
	Ro_SD1_SD2 = RO_SD1_SD2(Delay_Vth, Constraint_Vth);
	Vth_rand_sensitivity = Teta(VAR1, VAR2, Ro_SD1_SD2);//sqrt(pow(Constraiant_Vth_rand, 2) + pow(Delay_Vth_rand, 2) - (2 * ));
	VAR = Variance(Vth_sensitivity, Vth_rand_sensitivity);
	//******************
	
	Vth_rand_sensitivity = Constraiant_Vth_rand - Delay_Vth_rand;
	VAR = Variance(Vth_sensitivity, Vth_rand_sensitivity);

	//******************
	return (1 - Ta(Delay_mean, Constraint_mean, sqrt(Ro_SD1_SD2)));
	//return Ta((Constraint_mean - Delay_mean), sqrt(VAR));
	//return Ta(Constraint_mean, Delay_mean, sqrt(VAR));

	/*
	int i;
	double Vth_rand_sensitivity, VAR;
	vector<double> Vth_sensitivity;

	Vth_sensitivity.resize(16);
	for (i = 0; i < 16; i++)
		Vth_sensitivity[i] = Constraint_Vth[i] - Delay_Vth[i];
	Vth_rand_sensitivity = Constraiant_Vth_rand - Delay_Vth_rand;
	VAR = Variance(Vth_sensitivity, Vth_rand_sensitivity);
	//return Ta((Constraint_mean - Delay_mean) - (3 * sqrt(VAR)), 0, 1);
	
	//for (i = 0; i < 16; i++)
		//Vth_sensitivity[i] = Delay_Vth[i] - Constraint_Vth[i];
	//Vth_rand_sensitivity = Delay_Vth_rand - Constraiant_Vth_rand;
	//VAR = Variance(Vth_sensitivity, Vth_rand_sensitivity);
	//return Ta((Delay_mean - Constraint_mean) - (3 * sqrt(VAR)), 0, 1);
	
	return Ta(Constraint_mean, Delay_mean, sqrt(VAR));
	*/
}


double Yield_Computation(double Delay_mean, vector<double> &Delay_Vth, double Delay_Vth_rand, double Constraint)
{
	int i;
	double Vth_rand_sensitivity, VAR, VAR1, VAR2, Ro_SD1_SD2;
	vector<double> Vth_sensitivity;
	Engine *ep;

	//Normal_Distribution_CDF(Delay_mean / Variance(Delay_Vth, Delay_Vth_rand));
	if (!(ep = engOpen(NULL)))
	{
		cout << "error" << endl;
	}

	double xArray[1];
	xArray[0] = Constraint;
	mxArray* x = mxCreateDoubleMatrix(1, 1, mxREAL);
	memcpy((void *)mxGetPr(x), (void *)xArray, sizeof(xArray));
	engPutVariable(ep, "x", x);

	double muArray[1];
	muArray[0] = Delay_mean;
	mxArray* mu = mxCreateDoubleMatrix(1, 1, mxREAL);
	memcpy((void *)mxGetPr(mu), (void *)muArray, sizeof(muArray));
	engPutVariable(ep, "mu", mu);

	double sigmaArray[1];
	sigmaArray[0] = sqrt(Variance(Delay_Vth, Delay_Vth_rand));
	mxArray* sigma = mxCreateDoubleMatrix(2, 2, mxREAL);
	memcpy((void *)mxGetPr(sigma), (void *)sigmaArray, sizeof(sigmaArray));
	engPutVariable(ep, "sigma", sigma);

	engEvalString(ep, "p=normcdf(x,mu,sigma);");

	double *cresult;
	mxArray *mresult;
	mresult = engGetVariable(ep, "p");
	cresult = mxGetPr(mresult);

	return (*cresult);
}


void Yield_Optimization(vector<gate> &topological, vector<Critical> &Criticality, Sink *sink, int *count, double terminationCondition, double InitialCondition)
{
	int i = 0, j = 0, k;
	double temp = 0;
	double Constraint;
	*count = 0;
	//for (k = 0; k < 5; k++)
	while (true)
	{
		if (InitialCondition >= terminationCondition)
		{/*
			for (i = 0; i < 50; i++)
			{
				topological[Criticality[i].index].setWL_ratio(topological[Criticality[i].index].getWL_ratio() + 1);
				topological[Criticality[i].index].setDelay_mean();
				topological[Criticality[i].index].setAging_Delay_mean(10 * 365 * 24 * 60 * 60);
				topological[Criticality[i].index].setDelay_Vth_rand_sensitivity();
				topological[Criticality[i].index].setAging_Delay_Vth_rand_sensitivity(10 * 365 * 24 * 60 * 60);
				for (j = 0; j < 16; j++)
				{
					if (topological[Criticality[i].index].getDelay_Vth_sensitivity(j) != 0)
					{
						topological[Criticality[i].index].setDelay_Vth_sensitivity(j);
						topological[Criticality[i].index].setAging_Delay_Vth_sensitivity(j, 10 * 365 * 24 * 60 * 60);
					}
				}
			}

			(*count)++;
			SSTA(topological, sink, 0.98, 1.1, initial_optimization);
			temp = Variance(sink->Arrival_Vth_sensitivity, sink->Arrival_Vth_rand_sensitivity);
			cout << sink->Arrival_mean << "\t" << temp;
			Constraint = Yield_Computation(sink->Arrival_mean, sink->Arrival_Vth_sensitivity, sink->Arrival_Vth_rand_sensitivity, sink->initial_delay_mean,
				sink->initial_delay_Vth_sensitivity, sink->initial_delay_Vth_rand_sensitivity);
			cout << "\t" << Constraint << endl;
			*/

			break;
		}
		for (i = 0; (i <= 0.1 * Criticality.size()); i++)
		//for (i = 0; (Criticality[i].criticality >= 0.97 * Criticality[0].criticality); i++)
		//for (i = 0; i < 100; i++)
		{
			topological[Criticality[i].index].setWL_ratio(topological[Criticality[i].index].getWL_ratio() + 1);
			topological[Criticality[i].index].setDelay_mean();
			topological[Criticality[i].index].setAging_Delay_mean(10 * 365 * 24 * 60 * 60);
			topological[Criticality[i].index].setDelay_Vth_rand_sensitivity();
			topological[Criticality[i].index].setAging_Delay_Vth_rand_sensitivity(10 * 365 * 24 * 60 * 60);
			for (j = 0; j < 16; j++)
			{
				if (topological[Criticality[i].index].getDelay_Vth_sensitivity(j) != 0)
				{
					topological[Criticality[i].index].setDelay_Vth_sensitivity(j);
					topological[Criticality[i].index].setAging_Delay_Vth_sensitivity(j, 10 * 365 * 24 * 60 * 60);
				}
			}
		}

		(*count)++;
		SSTA(topological, sink, 0.6, 1.1, initial_optimization);
		temp = Variance(sink->Arrival_Vth_sensitivity, sink->Arrival_Vth_rand_sensitivity);
		cout << sink->Arrival_mean << "\t" << temp;
		Constraint = Yield_Computation(sink->Arrival_mean, sink->Arrival_Vth_sensitivity, sink->Arrival_Vth_rand_sensitivity, sink->RequiredArrivalTime);
		cout << "\t" << Constraint << endl;
		if (Constraint >= terminationCondition)
			break;
		Criticality_Computation(topological, sink);
		Rank_Criticality(topological, Criticality);
		temp = 0;
		i = 0;
	}
}


void GB_Gate_Sizing(vector<gate> &topological, vector<Age_Critical> &Aging_Criticality, vector<int> &inputs, Sink *sink, int *count, double terminationCondition, double Initial_Constraint)
{
	int i = 0, j = 0;
	double temp = 0;
	double Constraint;
	*count = 0;
	while (true)
	{
		if (Initial_Constraint <= terminationCondition)
			break;
		for (i = 0; (i <= (0.1 * Aging_Criticality.size())); i++)
		//for (i = 0; (Criticality[i].criticality >= 0.95 * Criticality[0].criticality); i++)
		//for (i = 0; i < 100; i++)
		{
			topological[Aging_Criticality[i].index].setWL_ratio(topological[Aging_Criticality[i].index].getWL_ratio() + 1);
			topological[Aging_Criticality[i].index].setDelay_mean();
			topological[Aging_Criticality[i].index].setAging_Delay_mean(10 * 365 * 24 * 60 * 60);
			topological[Aging_Criticality[i].index].setDelay_Vth_rand_sensitivity();
			topological[Aging_Criticality[i].index].setAging_Delay_Vth_rand_sensitivity(10 * 365 * 24 * 60 * 60);
			for (j = 0; j < 16; j++)
			{
				if (topological[Aging_Criticality[i].index].getDelay_Vth_sensitivity(j) != 0)
				{
					topological[Aging_Criticality[i].index].setDelay_Vth_sensitivity(j);
					topological[Aging_Criticality[i].index].setAging_Delay_Vth_sensitivity(j, 10 * 365 * 24 * 60 * 60);
				}
			}
		}

		(*count)++;
		SSTA(topological, sink, 0.6, 1.1, lifetime_optimization);
		temp = Variance(sink->Arrival_Vth_sensitivity, sink->Arrival_Vth_rand_sensitivity);
		cout << sink->Arrival_mean << "\t" << temp << endl;
		Constraint = sink->Arrival_mean + (3 * sqrt(temp));
		if (Constraint <= terminationCondition)
			break;
		Criticality_Computation(topological, sink);
		Aging_Computation(topological, sink);
		Pareto_Ranking(topological, sink, Aging_Criticality);
		temp = 0;
		i = 0;
	}
}


void GB_Gate_Sizing_Multiply(vector<gate> &topological, vector<Critical> &Aging_Criticality, vector<int> &inputs, Sink *sink, int *count, double terminationCondition, double Initial_Constraint)
{
	int i = 0, j = 0;
	double temp = 0;
	double Constraint;
	*count = 0;
	while (true)
	{
		if (Initial_Constraint <= terminationCondition)
			break;
		//for (i = 0; (i <= (0.1 * Aging_Criticality.size())); i++)
		for (i = 0; (Aging_Criticality[i].criticality >= 0.95 * Aging_Criticality[0].criticality); i++)
		//for (i = 0; i < 100; i++)
		{
			topological[Aging_Criticality[i].index].setWL_ratio(topological[Aging_Criticality[i].index].getWL_ratio() + 1);
			topological[Aging_Criticality[i].index].setDelay_mean();
			topological[Aging_Criticality[i].index].setAging_Delay_mean(10 * 365 * 24 * 60 * 60);
			topological[Aging_Criticality[i].index].setDelay_Vth_rand_sensitivity();
			topological[Aging_Criticality[i].index].setAging_Delay_Vth_rand_sensitivity(10 * 365 * 24 * 60 * 60);
			for (j = 0; j < 16; j++)
			{
				if (topological[Aging_Criticality[i].index].getDelay_Vth_sensitivity(j) != 0)
				{
					topological[Aging_Criticality[i].index].setDelay_Vth_sensitivity(j);
					topological[Aging_Criticality[i].index].setAging_Delay_Vth_sensitivity(j, 10 * 365 * 24 * 60 * 60);
				}
			}
		}

		(*count)++;
		SSTA(topological, sink, 0.6, 1.1, lifetime_optimization);
		temp = Variance(sink->Arrival_Vth_sensitivity, sink->Arrival_Vth_rand_sensitivity);
		cout << sink->Arrival_mean << "\t" << temp << endl;
		Constraint = sink->Arrival_mean + (3 * sqrt(temp));
		if (Constraint <= terminationCondition)
			break;
		Criticality_Computation(topological, sink);
		Aging_Computation(topological, sink);
		Pareto_Ranking_Multiply(topological, sink, Aging_Criticality);
		temp = 0;
		i = 0;
	}
}


void GB_Gate_Sizing_Critical_Aging(vector<gate> &topological, vector<Critical> &Aging_Criticality, vector<int> &inputs, Sink *sink, int *count, double terminationCondition, double Initial_Constraint)
{
	int i = 0, j = 0;
	double temp = 0;
	double Constraint;
	*count = 0;
	while (true)
	{
		if (Initial_Constraint <= terminationCondition)
			break;
		for (i = 0; (i <= (0.1 * Aging_Criticality.size())); i++)
		//for (i = 0; (Aging_Criticality[i].criticality >= 0.97 * Aging_Criticality[0].criticality); i++)
		//for (i = 0; i < 100; i++)
		{
			topological[Aging_Criticality[i].index].setWL_ratio(topological[Aging_Criticality[i].index].getWL_ratio() + 1);
			topological[Aging_Criticality[i].index].setDelay_mean();
			topological[Aging_Criticality[i].index].setAging_Delay_mean(10 * 365 * 24 * 60 * 60);
			topological[Aging_Criticality[i].index].setDelay_Vth_rand_sensitivity();
			topological[Aging_Criticality[i].index].setAging_Delay_Vth_rand_sensitivity(10 * 365 * 24 * 60 * 60);
			for (j = 0; j < 16; j++)
			{
				if (topological[Aging_Criticality[i].index].getDelay_Vth_sensitivity(j) != 0)
				{
					topological[Aging_Criticality[i].index].setDelay_Vth_sensitivity(j);
					topological[Aging_Criticality[i].index].setAging_Delay_Vth_sensitivity(j, 10 * 365 * 24 * 60 * 60);
				}
			}
		}

		(*count)++;
		SSTA(topological, sink, 0.6, 1.1, lifetime_optimization);
		temp = Variance(sink->Arrival_Vth_sensitivity, sink->Arrival_Vth_rand_sensitivity);
		cout << sink->Arrival_mean << "\t" << temp << endl;
		Constraint = sink->Arrival_mean + (2 * sqrt(temp));
		if (Constraint <= terminationCondition)
			break;
		Critical_Aging_Computation(topological, sink);
		Rank_Critical_Aging(topological, Aging_Criticality);
		temp = 0;
		i = 0;
	}
}


void GB_Gate_Sizing_Critical_Aging_itr(vector<gate> &topological, vector<Critical> &Aging_Criticality, vector<int> &inputs, Sink *sink, int *count, double terminationCondition, double Initial_Constraint)
{
	int itr, i = 0, j = 0;
	double temp = 0;
	double Constraint;
	*count = 0;
	for (itr = 0; itr < 11; itr++)
	{
		if (Initial_Constraint <= terminationCondition)
			break;
		//for (i = 0; (i <= (0.1 * Aging_Criticality.size())); i++)
		//for (i = 0; (Aging_Criticality[i].criticality >= 0.97 * Aging_Criticality[0].criticality); i++)
		for (i = 0; i < 100; i++)
		{
			topological[Aging_Criticality[i].index].setWL_ratio(topological[Aging_Criticality[i].index].getWL_ratio() + 1);
			topological[Aging_Criticality[i].index].setDelay_mean();
			topological[Aging_Criticality[i].index].setAging_Delay_mean(10 * 365 * 24 * 60 * 60);
			topological[Aging_Criticality[i].index].setDelay_Vth_rand_sensitivity();
			topological[Aging_Criticality[i].index].setAging_Delay_Vth_rand_sensitivity(10 * 365 * 24 * 60 * 60);
			for (j = 0; j < 16; j++)
			{
				if (topological[Aging_Criticality[i].index].getDelay_Vth_sensitivity(j) != 0)
				{
					topological[Aging_Criticality[i].index].setDelay_Vth_sensitivity(j);
					topological[Aging_Criticality[i].index].setAging_Delay_Vth_sensitivity(j, 10 * 365 * 24 * 60 * 60);
				}
			}
		}

		(*count)++;
		SSTA(topological, sink, 0.65, 1.1, lifetime_optimization);
		temp = Variance(sink->Arrival_Vth_sensitivity, sink->Arrival_Vth_rand_sensitivity);
		cout << sink->Arrival_mean << "\t" << temp << endl;
		Constraint = sink->Arrival_mean + (2 * sqrt(temp));
		//if (Constraint <= terminationCondition)
			//break;
		Critical_Aging_Computation(topological, sink);
		Rank_Critical_Aging(topological, Aging_Criticality);
		temp = 0;
		i = 0;
	}
}


void SUB_Delay_Distribution(double mean1, double var1, double mean2, double var2, double &mean_sub, double &var_sub)
{
	mean_sub = mean1 - mean2;
	var_sub = var1 - var2;
}