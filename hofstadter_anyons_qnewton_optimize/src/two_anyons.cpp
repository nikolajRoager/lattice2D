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
#include <chrono>
#include"minimize.hpp"

using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)


#include "approx.hpp"

int main(int argc, char* argv[])
{
    //Read arguments
    if (argc!=5)
    {
        cout<<"Invalid, argument, need "<<argv[0]<<" width height anyon_location_file potential_file"<<endl;
        return 1;
    }

    //If this fails, we get 0
    int nw = atoi(argv[1]);
    int nh = atoi(argv[2]);

    int N_sites = nh*nw;


    uint64_t N_particles = 4;

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

    if (N_sites %2 != 0 )
    {
        cout <<"Need even number of sites"<<endl;
        return 1;
    }

    if (0 > N_particles)
    {
        cout<<"Negative particle number "<< N_particles<<" not allowed"<<endl;
        return 1;
    }


    vector<uint64_t> anyon_sites(N_sites);//The sites containing 0: no anyons, 1: first anyon 2: second anyon

    ifstream anyon_file(argv[3]);

    if (!anyon_file.is_open())
    {
        cout<<"Could not open "<<argv[3]<<endl;
        return 1;
    }

    for (uint64_t i = 0; i < N_sites; ++i)
    {
        uint64_t I =0;
        if (! (anyon_file>>I))
        {
            cout<<argv[3]<<" did not match number of sites (site "<<i<<" could not be read)" <<endl;
            return 1;
        }

        if (I>2)
            cout<<"Site "<<i<<" could not be read in "<<argv[3]<<" is outside range, should be 0 (no anyons), 1 or 2 (either anyon), got "<<I<<endl;

        anyon_sites[i]= I;
    }


    anyon_file.close();


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

    vec lattice_dens0= get_density(states0,prop0_list, N_states0, nw, nh);
    //print_density_data(lattice_dens0,  nw, nh);//Ok now print it so gnuplot can plot it
    //cerr << endl;


    //Now add in two anyons, removing one particle (so keeping phi the same) and adding a potential

    //But first. Set up all the variables we would like to save, and need to use

    //Traping potentials AND auxillary potential on all other sites put together
    uint64_t half_N_sites=N_sites/2;
    vec V(half_N_sites);


    ifstream potential_file(argv[4]);

    if (!potential_file.is_open())
    {
        cout<<"Could not open "<<argv[4]<<endl;
        return 1;
    }

    for (uint64_t i = 0; i < half_N_sites; ++i)
    {
        if (! (potential_file>>V[i]))
        {
            cout<<"Potential in "<<argv[4]<<" did not match number of sites (potential at site "<<i<<" could not be read)"<<endl;
            return 1;
        }
    }


    potential_file.close();

    vec lattice_dens1;//Will be written to, if we want to plot the actual density

    //The actual density difference we found
    vec q(N_sites);


    //I will use a lambda function, to make getting the fitness, densities and charge difference easier, though C++ allows me to modify the captures I NEVER DO THAT IN A PROJECT WHERE I ALSO USE MULTITHREADING, FOR OBVIOUS REASONS, the things actually written to are all given as arguments
    auto get_fitness =[&lattice_dens0,&anyon_sites,N_states1, &states1,N_sites,N_particles,nw,nh](vec V, vec& lattice_dens1, vec& q) -> double
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
            states1,    //written to
            a2_eigval,  //written to
            a2_eigvec,   //written to
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


        //Check

        //
        q=vec(N_sites);
        double fitness=0;

        double anyon0_charge=0;
        double anyon1_charge=0;

        for (uint64_t i = 0; i < N_sites; ++i)
        {
                double n0 =lattice_dens0(i);
                double n1 =lattice_dens1(i);
                q(i)=n0-n1;
                if (anyon_sites[i]==0)
                    fitness+=std::abs(q(i));//Whatever charge is here IS WRONG
                else if (anyon_sites[i]==1)//But I don't care where the charge of the anyons are, as long as thje charge in total adds up to 0.5
                    anyon0_charge+=std::abs(q(i));
                else if (anyon_sites[i]==2)
                    anyon1_charge+=std::abs(q(i));
        }

        fitness+=std::abs(anyon0_charge-0.5)+std::abs(anyon1_charge-0.5);



        return fitness;
    };

    auto get_fitness_sym =[&get_fitness,half_N_sites,N_sites](vec _V, vec& lattice_dens1, vec& q) -> double
    {
        vec V(N_sites);
        for (uint64_t i = 0; i < half_N_sites; ++i)
        {
            V(i) = _V(i);
            V(N_sites-i-1)=_V(i);
        }

        return get_fitness(V,lattice_dens1,q);
    };

    cout<<"Measuring time of get fitness function, by running it 50 times "<<endl;
    double time=0;
    double fitness;
    auto t1 = std::chrono::high_resolution_clock::now();
    for (uint64_t i = 0; i < 50; ++i)
    {
        fitness = get_fitness_sym (V,lattice_dens1,q);

    }
    auto t2 = std::chrono::high_resolution_clock::now();

    /* Getting number of milliseconds as a double. */
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;
    time += ms_double.count();

    cout<<"Got "<<time/50<<" ms per iteration "<<endl;

    print_density_data(q , nw, nh);
    cerr << endl;

    cout<<"Fitness with unoptimized potential "<<argv[4]<<" : "<<fitness<<endl;


    size_t steps;
    cout<<std::setprecision(16);

    auto t_opt_start = std::chrono::high_resolution_clock::now();

    //Optimize symmetric potential
    vec V_optimized_half = qnewton([&get_fitness_sym](vec V)
    {
        //Don't output these things, but create them locally so that we don't run into problems with multithreading
        vec lattice_dens1;
        vec q;
        return get_fitness_sym(V,lattice_dens1,q);
    },V,steps,128,1e-5,true);




    auto t_opt_end = std::chrono::high_resolution_clock::now();
    std::chrono::minutes duration = std::chrono::duration_cast<std::chrono::minutes>(t_opt_end - t_opt_start);
    cout<<"Got optimized symmetric potential in "<<steps<<" steps ("<<duration.count()<<" minutes) with "<<get_fitness_sym(V_optimized_half,lattice_dens1,q)<<endl;

    vec V_optimized(N_sites);
    for (uint64_t i = 0; i < half_N_sites; ++i)
    {
        V_optimized(i) = V_optimized_half(i);
        V_optimized(N_sites-i-1)=V_optimized_half(i);
    }
/*
    cout<<"Now turning off symmetry to optimize a little better "<<endl;

    //This is much slower
    t_opt_start = std::chrono::high_resolution_clock::now();
    V_optimized = qnewton([&get_fitness](vec V)
    {
        //Don't output these things, but create them locally so that we don't run into problems with multithreading
        vec lattice_dens1;
        vec q;
        return get_fitness(V,lattice_dens1,q);
    },V_optimized,steps,16,1e-5,true);
    t_opt_end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::minutes>(t_opt_end - t_opt_start);

    cout<<"Got optimized potential in another "<<steps<<" steps ("<<duration.count()<<" minutes) with "<<get_fitness(V_optimized,lattice_dens1,q)<<endl;
*/
    //q is written to as soon as we call get_fitness
    print_density_data(q , nw, nh);
    cerr << endl;


    print_density_data(V_optimized, nw, nh);
    cerr << endl;

    cout<<endl<<"potential:"<<endl;
    cout<<setprecision(16);
    for (double Vi : V_optimized)
        cout<<Vi<<'\n';
    cout<<endl;
    return 0;
}
