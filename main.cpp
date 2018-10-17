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

double energy(vector<vector<int>> &v, int rows, int cols)
{
	double energy_value = 0;
	double j_const = 1;

	/* Convert the 0's in the matrix to -1's */
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			if(v[i][j] == 0)
				v[i][j] = -1;
		}
	}
	
	

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
	if(v[spin_y][spin_x] == -1)
	{
		v[spin_y][spin_x] = 1;
	}
	else
	{
		v[spin_y][spin_x] = -1;
	}
}

int main()
{
	
	string input = "01100101 01100001 10010110 11110000 00001010 11110101 10111001 11000001";
	cout << "The initialized state is " << input << endl;
	vector<vector<int>> spins, new_spins;

	cout << "Its matrix is " << endl;


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

	int comparison = dis(gen);

	/* Convert the input string into a matrix */
	spins.resize(cols);
	spins = from_string(input,spins);
	print_state(spins);

	cout << endl;

	new_spins = spins;
	change_spins(new_spins,spin_x,spin_y);
	print_state(new_spins);
	
	double beta = 1.5;

	int flag = 0;

	

	cout << endl;
	cout << "Writing the new state..." << endl;

	ofstream out_file;
	out_file.open("out.txt");

	for(int sweep = 0; sweep < 10000; sweep++)
	{
		for(int i = 0; i < rows*cols; i++)
		{
			spin_x = distr(generator);
			spin_y = distr(generator);

				double energy_value_old =  energy(spins, rows, cols);
				double energy_value_new =  energy(new_spins, rows, cols);
				double delta_energy_new_old = delta_e(spins,rows,cols,0,0,1);
				double delta_energy_old_new = delta_e(spins,rows,cols,0,0,0);



			if((delta_energy_new_old/delta_energy_old_new)*(exp(-beta*energy_value_new)/exp(-beta*energy_value_old)) > comparison)
			{
				change_spins(spins,spin_x,spin_y);
			}
		}

		flag++;

		if(flag >= 10)
		{
			for(size_t i = 0; i < spins.size(); i++)
			{
					for(size_t j = 0; j < spins[i].size(); j++)
					{
						if(spins[i][j] == -1)
							spins[i][j] = 0;

						out_file << spins[i][j] << " ";
					}
					out_file << endl;
			}	
		}	
	}
	out_file.close();
/*
	cout << "The randomly selected spin to switch is located at ";
	cout << "[" << spin_x << "," << spin_y << "]" << endl;

	cout << "The energy of the current state is " << energy_value_old << endl;
	cout << "The difference of energy between the two states is " << delta_energy_new_old << endl;

*/
	cout << "The new state is " << endl;
	print_state(spins);

	return 0; 
}