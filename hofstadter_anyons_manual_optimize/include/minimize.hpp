#pragma once
#include<cstdint>
#include<iostream>
#include<string>
#include<armadillo>
#include<functional>

using state = uint64_t;
using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)


vec qnewton( function<double(vec)> func , vec x0, size_t& steps, double acc=1e-2, bool verbose=false);

//Mirror function to maximize instead
vec max_qnewton( function<double(vec)> func , vec x0, size_t& steps, double acc=1e-2, bool verbose=false);

//Versions without counting number of steps
vec qnewton( function<double(vec)> func , vec x0, double acc=1e-2, bool verbose=false);

//Mirror function to maximize instead
vec max_qnewton( function<double(vec)> func , vec x0, double acc=1e-2, bool verbose=false);
