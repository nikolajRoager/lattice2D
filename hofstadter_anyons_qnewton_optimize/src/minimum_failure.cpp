#include<iostream>
#include<complex>
#include<exception>
#include <armadillo>
#include<thread>
#include<string>

using namespace std;
using namespace arma;

#include "approx.hpp"

std::mutex stupid_lock;



void get_eigenvalue(int id, size_t N, double* out)
{
    sp_cx_mat H (N, N);//Sparse matrix, as we have mostly 0 entries, initialized as only zeros by default

    //As a demonstration, a hermitian matrix, which should have real eigenvalues
    for (uint i = 0; i <N; ++i)
    {
        H(i, i)=1*i;

        if (i>0)
            H(i, i-1)=-1i*double(i);
        if (i<N-1)
            H(i, i+1)=+1i*double(i);
    }


    //Somewhere to put the eigenvalues and vectors
    cx_vec eigval;
    cx_mat eigvec;

    //stupid_lock.lock();
    //It makes no sense that I need thism everything was created locally in this function
    eigs_gen( eigval, eigvec, H, 1, "sr");//The lowest eigenvalue
    //stupid_lock.unlock();

    *out = eigval[0].real();//Should be real, as the matrix is hermitian

    if ( eigval[0].imag()> 1e-9 )//Sanity check, if this is not real something has failed (in practice allow for slight rounding errors)
    {
        throw std::runtime_error("ERROR in thread "+to_string(id)+", energy "+to_string(eigval[0].real())+"+"+to_string(eigval[0].imag())+"*i + has non-zero imaginary part");
    }
}

int main()
{
    size_t N = 20000;//The error only happens if the functions take a substantial amount of time, this takes around half a second to diagonalise on my computer
    double out1=0;
    double out2=0;


    thread t1(get_eigenvalue,1,N,&out1);
    thread t2(get_eigenvalue,2,N,&out2);

    t1.join();
    t2.join();
    cout<<"These numbers should be the same "<<out1<<" "<<out2<<endl;

    return 0;
}
