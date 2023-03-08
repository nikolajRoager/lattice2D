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
#include<fstream>

//The actual math library
#include <armadillo>
#include"minimize.hpp"

using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)


#include "approx.hpp"

int main(int argc, char* argv[])
{
    //Read arguments
    if (argc!=8)
    {
        cout<<"Invalid, argument, need "<<argv[0]<<" width height x0 y0 x1 y1 potential_file"<<endl;
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


    double empty_eigval;
    cx_vec empty_eigvec;

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

    //To calculate the density of particles in the lattice, we need the probability of being in each basis state
    vector<double> prop0_list(N_states0);
    for (uint64_t j = 0; j < N_states0; ++j)
    {
        prop0_list[j]= norm(empty_eigvec[j]);
    }

    mat lattice_dens0= get_density(states0,prop0_list, N_states0, nw, nh);
    //print_density_data(lattice_dens0,  nw, nh);//Ok now print it so gnuplot can plot it
    //cerr << endl;


    //Now add in two anyons, removing one particle (so keeping phi the same) and adding a potential

    //But first. Set up all the variables we would like to save, and need to use

    //Traping potentials AND auxillary potential on all other sites put together
    vec V(N_sites);


    ifstream potential_file(argv[7]);

    if (!potential_file.is_open())
    {
        cout<<"Could not open file with potentials"<<endl;
        return 1;
    }

    for (uint64_t i = 0; i < N_sites; ++i)
    {
        if (! (potential_file>>V[i]))
        {
            cout<<"Potential file did not match number of sites"<<endl;
            return 1;
        }
    }


    potential_file.close();

    mat lattice_dens1;//Will be written to, if we want to plot the actual density

    //This is what we really want, two perfectly localized anyons
    mat q_optimum=mat(nw,nh);
    q_optimum(x0,y0)=0.5;
    q_optimum(x1,y1)=0.5;


    //Tje actual density difference we found
    mat q(nw,nh);

    //I will use a lambda function, to make getting the fitness, densities and charge difference easier
    auto get_fitness =[&lattice_dens0,&lattice_dens1 ,&q_optimum,&q,&states1,N_states1,N_sites,N_particles,nw,nh](vec V) -> double
    {

        double a2_eigval;
        cx_vec a2_eigvec;

        get_state(
            N_states1,
            N_sites,
            N_particles-1,
            2,//N_anyons,
            nw,
            nh,
            states1,
            a2_eigval,
            a2_eigvec,
            V
        );
        //Get the propability of finding the system in either basis state, this is one step before finding the expected number of particles in each site

        vector<double> prop1_list(N_states1);
        //I strictly don't need to use uint64 here, but the compiler just gets real mad if I use different types
        for (uint64_t j = 0; j < N_states1; ++j)
        {
            prop1_list[j]= norm(a2_eigvec[j]);
        }

        lattice_dens1 = get_density(states1,prop1_list, N_states1, nw, nh);

        //
        q=mat(nw,nh);
        double fitness=0;
        for (uint64_t y = 0; y<nw; ++y)
        {
            for (uint64_t x = 0; x<nw; ++x)
            {
                double n0 =lattice_dens0(x,y);
                double n1 =lattice_dens1(x,y);
                q(x,y)=n0-n1;
                fitness += std::abs((n0-n1)-q_optimum(x,y));

            }
        }

        return fitness;
    };




    double fitness = get_fitness (V);

    print_density_data(q , nw, nh);
    cerr << endl;

    cout<<"Fitness with unoptimized potential "<<argv[7]<<" : "<<fitness<<endl;


    size_t steps;
    cout<<std::setprecision(16);
    vec V_optimized = qnewton(get_fitness,V,steps,1e-5,true);

    cout<<"Got optimized potential in "<<steps<<" with "<<get_fitness(V_optimized)<<endl;

    //q is written two as soon as we call get_fitness
    print_density_data(q , nw, nh);
    cerr << endl;



    return 0;
}
