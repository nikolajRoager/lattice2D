#pragma once

//Has already been done in my header, but I like to include everything explicitly in every file (the compiler will ignore it)
#include<cstdint>

//NOTE! std::vector is a dynamic sized list, it is not a vector in the linear algebra sense of the word, that will be the class vec from the armadillo library
#include<vector>
#include<complex>

//The actual math library
#include <armadillo>

using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)

void get_state(
uint64_t N_states,
uint64_t N_sites,
uint64_t N_particles,
uint64_t N_anyons,
uint64_t nw,
uint64_t nh,
const vector<state_t>& states,
double& this_eigval,
cx_vec& this_eigvec,
vec& V);

//Default version without a potential
void get_state(
uint64_t N_states,
uint64_t N_sites,
uint64_t N_particles,
uint64_t N_anyons,
uint64_t nw,
uint64_t nh,
const vector<state_t>& states,
double& this_eigval,
cx_vec& this_eigvec
);
