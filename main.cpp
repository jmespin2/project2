#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>

#include "input.h"


using namespace std;

/*
InputClass read_file;
read_file.Read("input_file.txt");

double beta = read_file.toDouble(read_file.GetVariable("beta"));
int Lx = read_file.toInteger(read_file.GetVariable("Lx"));
int Ly = read_file.toInteger(read_file.GetVariable("Ly"));
*/

void print_state(vector<vector<int>> &state)
{

	for(size_t i = 0; i < state.size(); i++)
	{
		for(size_t j = 0; j < state[i].size(); j++)
		{
			cout << state[i][j] << " ";
		}
		cout << endl;
	}

	return;
}
/*
void write_state(vector<vector<int>> &state)
{

	for(size_t i = 0; i < state.size(); i++)
	{
		for(size_t j = 0; j < state[i].size(); j++)
		{
			if(state[i][j] == -1)
				state[i][j] = 0;

			out_file << state[i][j] << " ";
		}
		out_file << endl;
	}
	out_file.close();

	return;
}
*/
vector<vector<int>> from_string(string &input, vector<vector<int>> &v)
{
	int j = 0;

	for(size_t i = 0; i < input.length(); i++)
	{
		if(input[i] == ' ')
		{
			j++;
		}

		else
		{
			v[j].push_back(input[i] - 48);					//convert ASCII to int
		}

		
	}


	return v;
}

void convert_0s_to_1s(vector<vector<int>> &v, int rows, int cols)
{
		/* Convert the 0's in the matrix to -1's */
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			if(v[i][j] == 0)
				v[i][j] = -1;
		}
	}
}

double energy(vector<vector<int>> &v, int rows, int cols)
{
	double energy_value = 0;
	double j_const = 1;

	convert_0s_to_1s(v,rows,cols);	

	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			int current_spin = v[i][j];
			
			/* Check the up element's boundary */
			if((i - 1) < 0)
			{
				energy_value = energy_value + (j_const*current_spin*v[rows-1][j]);
			}
			else
			{
				energy_value = energy_value + (j_const*current_spin*v[i-1][j]);	
			}

			/* Check the right element's boundary */
			if((j + 1) > (cols - 1))
			{
				energy_value = energy_value + (j_const*current_spin*v[i][0]);
			}
			else
			{
				energy_value = energy_value + (j_const*current_spin*v[i][j+1]);
			}
		}

	}

	energy_value = -energy_value;

	return energy_value;
}

double delta_e(vector<vector<int>> &v, int rows, int cols, int spin_x, int spin_y, int arg)
{
	double energy_value = energy(v,rows,cols);
	double old_energy_value = energy_value;
	double j_const  = 1;

	convert_0s_to_1s(v,rows,cols);

	/* Up */
	if((spin_y - 1) < 0)
	{
		energy_value = energy_value + (j_const*(-v[spin_y][spin_x])*v[rows-1][spin_x]);
		energy_value = energy_value - (j_const*(v[spin_y][spin_x])*v[rows-1][spin_x]);
	}
	else
	{
		energy_value = energy_value + (j_const*(-v[spin_y][spin_x])*v[spin_y-1][spin_x]);
		energy_value = energy_value - (j_const*v[spin_y][spin_x]*v[spin_y-1][spin_x]);
	}

	/* Right */
	if((spin_x + 1) > (cols - 1))
	{
		energy_value = energy_value + (j_const*(-v[spin_y][spin_x])*v[spin_y][0]);
		energy_value = energy_value - (j_const*(v[spin_y][spin_x])*v[spin_y][0]);
	}
	else
	{
		energy_value = energy_value + (j_const*(-v[spin_y][spin_x])*v[spin_y][spin_x+1]);
		energy_value = energy_value - (j_const*v[spin_y][spin_x]*v[spin_y][spin_x+1]);
	}

	/* Down */
	if((spin_y + 1) > (rows - 1))
	{
		energy_value = energy_value + (j_const*(-v[spin_y][spin_x])*v[0][spin_x]);
		energy_value = energy_value - (j_const*(v[spin_y][spin_x])*v[0][spin_x]);
	}
	else
	{
		energy_value = energy_value + (j_const*(-v[spin_y][spin_x])*v[spin_y+1][spin_x]);
		energy_value = energy_value - (j_const*v[spin_y][spin_x]*v[spin_y+1][spin_x]);
	}

	/* Left */
	if((spin_x - 1) < 0)
	{
		energy_value = energy_value + (j_const*(-v[spin_y][spin_x])*v[spin_y][cols-1]);
		energy_value = energy_value - (j_const*(v[spin_y][spin_x])*v[spin_y][cols-1]);
	}
	else
	{
		energy_value = energy_value + (j_const*(-v[spin_y][spin_x])*v[spin_y][spin_x]);
		energy_value = energy_value - (j_const*v[spin_y][spin_x]*v[spin_y][spin_x]);
	}

	if(arg == 0)
		return old_energy_value - energy_value;
	else
		return energy_value - old_energy_value;
}


void change_spins(vector<vector<int>> &v, int spin_x, int spin_y)
{
	/* Change the spin. */
	if(v[spin_y][spin_x] == 0)
	{
		v[spin_y][spin_x] = 1;
	}
	else if(v[spin_y][spin_x] == -1)
	{
		v[spin_y][spin_x] = 1;
	}
	else
	{
		v[spin_y][spin_x] = 0;
	}

}

void find_probability_sum(double &P, double energy, double beta)
{
	double probability = exp(-beta*energy);
	//cout << "probability is " << probability << endl;
	//cout << "energy is " << energy << endl;
	P += probability;
}

void write_probabilities(double &P, vector<double>E)
{
	ofstream probability_file;
	probability_file.open("probability.txt");

	for(size_t i = 0; i < E.size(); i++)
	{
		probability_file << abs(E[i]/P) << endl;
	}
	//cout << P << endl;
	probability_file.close();
}

double find_magnetization(vector<vector<int>> &v, int sweeps)
{
	double M = 0;

	return M;
}

int main()
{
	

	
	double beta = 0.3;
	/* 27x27 input */
	//string input = "100010110111010110101010110 100010110100010110101010110 101110110111010110101010110 100110110100010110101010110 100010110111011111101010110 101010110111010110101010110 000010110111010111111010111 100010110111011110101010110 111110110111010110101011111 101110110111010110101010110 100010110111010110101010110 011010110111010110101110110 100010100001010110101010110 100010110100010110101010110 100010000111010110101010110 100010110110010110101010110 100010110111010110111110110 101110110111010110101010110 111110110111010111101110110 101010110100010110101010000 100010110000010110101010110 100010110111010110101010000 000010110001010010101010110 100010110100010110101010110 000010110110010100101010111 100011110111011110101000110 100010111111011110111010101";
	
	string input = "111 111 111";

	//cout << "The initialized state is " << input << endl;
	vector<vector<int>> spins, new_spins;
	vector<double> saved_energies;
	double P = 0;
	double E = 0;


	/* Check the number of spaces */
	int cols = 0;
	int rows = 1;

	for(size_t i = 0; i < input.length(); i++)
	{
		if(input[i] == ' ')
		{
			rows++;
		}

	}

	for(size_t i = 0; i < input.length(); i++)
	{
		
		if(input[i] == ' ')
		{
			break;
		}
		cols++;

	}


	/* Random number generator here */
	const int range_from = 0;
	const int range_to = rows - 1;
	const int range_to_2 = 1;

	random_device rand_dev;
	mt19937 generator(rand_dev());
	uniform_int_distribution<int> distr(range_from, range_to);

	random_device rand_dev2;
	mt19937 gen(rand_dev2());
	uniform_int_distribution<int> dis(range_from, range_to_2);

	int spin_x = distr(generator);
	int spin_y = distr(generator);

	
	/* Convert the input string into a matrix */
	spins.resize(cols);
	spins = from_string(input,spins);
	
	cout << rows << ' ' << cols << endl;
	print_state(spins);

	cout << endl;

	new_spins = spins;
	
	int flag = 0;

	cout << endl;
	cout << "Writing the new state..." << endl;
	ofstream out_file;
	out_file.open("out.txt");


	int probability_flag = 0;
	/* Sweep here */
	for(int sweep = 0; sweep < 10000; sweep++)
	{
		for(int i = 0; i < rows*cols; i++)
		{
			spin_x = distr(generator);
			spin_y = distr(generator);
			//cout << spin_x << " " << spin_y << endl;

			/* Create a temporary configuration to check the energy */
			change_spins(new_spins,spin_x,spin_y);

			//double energy_value_old =  energy(spins, rows, cols);
			//double energy_value_new =  energy(new_spins, rows, cols);
			//double delta_energy_old_new = delta_e(spins,rows,cols,spin_x,spin_y,1);
			double delta_energy_new_old = delta_e(spins,rows,cols,spin_x,spin_y,0);

			double e_compare = exp(-beta*delta_energy_new_old);
			int comparison = dis(gen);

			/* Decide if we accept the new configuration */
			if(e_compare > comparison)
			{
				spins = new_spins;				
			}

			/* Calculate per configuration */
			saved_energies.push_back(exp(-beta*E));
			find_probability_sum(P,E,beta);
	
		}


			/* Take a snapshot of the sweep */
			if(flag >= 10)
			{
				for(size_t k = 0; k < spins.size(); k++)
				{
						for(size_t j = 0; j < spins[k].size(); j++)
						{
							if(spins[k][j] == -1)
								spins[k][j] = 0;

							out_file << spins[k][j];
							E = energy(spins,rows,cols);
						}
						out_file << endl;						
				}	
				
			}
			flag++;
			
			//cout << "P is " << P << endl;
	}

	write_probabilities(P,saved_energies);

	
	out_file.close();
/*
	cout << "The randomly selected spin to switch is located at ";
	cout << "[" << spin_x << "," << spin_y << "]" << endl;

	cout << "The energy of the current state is " << energy_value_old << endl;
	cout << "The difference of energy between the two states is " << delta_energy_new_old << endl;

*/

	/* Convert the 0's in the matrix to -1's */
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			if(spins[i][j] == -1)
				spins[i][j] = 0;
		}
	}
	

	cout << "The new state is " << endl;
	print_state(spins);

	return 0; 
}