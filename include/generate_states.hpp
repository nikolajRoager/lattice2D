#pragma once
//The function for generating the states, and some of the math required to do that

#include<cstdint>
#include<vector>

//It is critically important that I use the same integer type for all binary operations, and seeing that the state is 64 bit, then almost everything else need to be 64 bit as well.
using state = uint64_t;

using namespace std;

uint64_t binSearch(const vector<state>& states, state target, uint64_t low=0 ,uint64_t high=-1 /*Intentional underflow, is truncated to the end of the list automatically */);

uint64_t bin_coef(uint64_t n, uint64_t m);


//Print this state (human readable), just as a sanity check
void print_state(state S,  int w, int h);//Assumes little Endian states

//Print this state in a format gnuplot can plot
void print_state_data(state S,  int w, int h);//Assumes little Endian states

//Print this state in a format gnuplot can plot
void print_superpositions_data(const vector<state>& states,vector<double> superpositions , uint64_t& n_states,  int w, int h);//Assumes little Endian states


//Generate all states and stores them here, assumes that size has been allocated correctly
void generate(vector<state>& states,
              uint64_t& N_states,
              uint64_t N_sites,
              uint64_t N_particles);
