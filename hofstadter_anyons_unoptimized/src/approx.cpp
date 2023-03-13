#include "approx.hpp"

#include<vector>
#include<iostream>
#include<string>
#include<complex>
#include<exception>

//The actual math library
#include <armadillo>



using namespace std;

bool approx(double a, double b,double tau, double epsilon)
{

    if (abs(a-b)<=tau)
        return true;
    else if (abs(a-b)/(abs(a)+abs(b))<=epsilon)
        return true;
    else
        return false;
}

