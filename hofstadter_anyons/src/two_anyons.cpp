#include"generate_states.hpp"
#include"get_state.hpp"

//Has already been done in my header, but I like to include everything explicitly in every file (the compiler will ignore it)
#include<cstdint>

//NOTE! std::vector is a dynamic sized list, it is not a vector in the linear algebra sense of the word, that will be the class vec from the armadillo library
#include<vector>
#include<iostream>
#include<string>
#include<complex>
#include<exception>
#include<iomanip>

//The actual math library
#include <armadillo>

using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)


#include "approx.hpp"

int main(int argc, char* argv[])
{
    //Read arguments
    if (argc!=7)
    {
        cout<<"Invalid, argument, need "<<argv[0]<<" width height x0 y0 x1 y1"<<endl;
        return 1;
    }

    //If this fails, we get 0
    int nw = atoi(argv[1]);
    int nh = atoi(argv[2]);

    int N_sites = nh*nw;

    int x0 = atoi(argv[3]);
    int y0 = atoi(argv[4]);
    int x1 = atoi(argv[5]);
    int y1 = atoi(argv[6]);


    uint64_t N_particles = 3;

    if (nw==0)
    {
        cout << "Invalid argument "<<argv[1]<<", must be integer"<<endl;
        return 1;
    }
    if (nh==0)
    {
        cout << "Invalid argument "<<argv[2]<<", must be integer"<<endl;
        return 1;
    }
    if (N_sites >sizeof(state_t)*8)
    {
        cout <<"This program only allows up to "<<(sizeof(state_t)*8)<<" sites, got  "<<nw<<"x"<<nh<<"="<<N_sites<<endl;
        return 1;
    }

    if (N_sites < N_particles)
    {
        cout <<"More particles "<<N_particles<<" than lattice sites  "<<nw<<"x"<<nh<<"="<<N_sites<<endl;
        return 1;
    }


    if (0 > N_particles)
    {
        cout<<"Negative particle number "<< N_particles<<" not allowed"<<endl;
        return 1;
    }


    cout << nw<<"x"<<nh<<" Lattice with "<<N_particles<<" particles for "<<N_sites<<" sites"<<endl;

    cout << "Generating and printing ALL legal states to stderr"<<endl;



    if (x0 >= nw)
    {
        cout<<"anyon x location "<< x0<<" outside lattice, bust be from 0 to "<<(nw-1)<<endl;
        return 1;
    }

    if (y0 >= nw)
    {
        cout<<"anyon x location "<< x0<<" outside lattice, bust be from 0 to "<<(nw-1)<<endl;
        return 1;
    }

    if (x1 >= nw)
    {
        cout<<"anyon x location "<< x0<<" outside lattice, bust be from 0 to "<<(nw-1)<<endl;
        return 1;
    }
    if (y1 >= nw)
    {
        cout<<"anyon x location "<< x0<<" outside lattice, bust be from 0 to "<<(nw-1)<<endl;
        return 1;
    }


    //I am very unlikely to get a number of states which needs uint64_t here, but I do need to use the same type as I use for the states (written as a 64 bit number), some of the binary operations may fail otherwise

    //Basis states with 3 or 2 particles
    uint64_t N_states0=0;
    vector<state_t> states0(N_states0);

    uint64_t N_states1=0;
    vector<state_t> states1(N_states1);

    generate(states0,N_states0,N_sites,N_particles);
    generate(states1,N_states1,N_sites,N_particles-1);

    cout <<"Found "<<N_states0<<" states with "<<N_particles<<" particles"<<endl;
    cout <<"And "<<N_states1<<" states with "<<(N_particles-1)<<" particles"<<endl;

    double empty_eigval;
    cx_vec empty_eigvec;

    cout <<"Getting empty ground state"<<endl;
    get_state(
    N_states0,
    N_sites,
    N_particles,
    0,//N_anyons,
    nw,
    nh,
    states0,
    empty_eigval,
    empty_eigvec);

    //Print this state density to a table, which gnuplot can plot
    vector<double> prop0_list(N_states0);
    for (uint64_t j = 0; j < N_states0; ++j)
    {
        prop0_list[j]= norm(empty_eigvec[j]);
    }

    print_superpositions_data(states0,prop0_list , N_states0,  nw, nh);
    cerr << endl;


    //Now add in two anyons, removing one particle (so keeping phi the same) and adding a potential


    vec V(N_sites);

    V[x0+y0*nw]=1;
    V[x1+y1*nw]=1;
    cout <<"Getting ground state with potential and 1 fewer particles"<<endl;

    get_state(
    N_states1,
    N_sites,
    N_particles-1,
    2,//N_anyons,
    nw,
    nh,
    states1,
    empty_eigval,
    empty_eigvec,
    V);

    //Get the propability of finding the system in either basis state, this is one step before finding the expected number of particles in each site

    vector<double> prop1_list(N_states1);
    //I strictly don't need to use uint64 here, but the compiler just gets real mad if I use different types
    for (uint64_t j = 0; j < N_states1; ++j)
    {
        prop1_list[j]= norm(empty_eigvec[j]);
       // prop2_list[j]= prop0_list[j]-prop1_list[j];
    }

    print_superpositions_data(states1,prop1_list , N_states1,  nw, nh);
    cerr << endl;

    //mat lattice_charge = get_density(states,prop2_list, N_states, nw, nh);
    mat lattice_dens0= get_density(states0,prop0_list, N_states0, nw, nh);
    mat lattice_dens1 = get_density(states1,prop1_list, N_states1, nw, nh);


    std::cout << std::fixed;
    std::cout << std::setprecision(3);
    double sum0 = 0;
    double sum1 = 0;
    double sum2 = 0;
    cout<<"Density at each site:"<<endl;
    cout<<" |\t";
    for (uint64_t x = 0; x<nw; ++x)
        cout<<x<<"\t";
    cout<<'\n';

    mat q(nw,nh);
    for (uint64_t y = 0; y<nw; ++y)
    {
        cout<<y<<"|\t";
        for (uint64_t x = 0; x<nw; ++x)
        {
            double n0 =lattice_dens0(x,y);
            double n1 =lattice_dens1(x,y);
            sum0 += n0;
            sum1 += n1;
            q(x,y)=n0-n1;
            sum2 += q(x,y);

            cout<<q(x,y)<<"\t";
        }
        cout<<'\n';
    }
    cout<<flush;


    print_density_data(q , nw, nh);
    cerr << endl;


    cout<<"sum density 0 = "<<sum0<<endl;
    cout<<"sum density 1 = "<<sum1<<endl;
    cout<<"sum charge = "<<sum2<<endl;
    return 0;
}
