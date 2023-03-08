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

//The actual math library
#include <armadillo>

using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)


#include "approx.hpp"


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
)
{
vec V(N_sites);
get_state(
N_states,
N_sites,
N_particles,
N_anyons,
nw,
nh,
states,
this_eigval,
this_eigvec,
V);
}
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
vec& V)
{
    //Make our hamiltonian matrix
    sp_cx_mat H (N_states, N_states);//Sparse matrix, as we have mostly 0 entries, initialized as only zeros by default

    //Loop through all our states and calculate all matrix elements
    //In my case, I know only nearest neighbours interact (Hofstadter model), so I only loop through each state once, and I only consider neighbour with a lower state than this state (Lower triangle, we know the matrix is hermitian! so we will just mirror everything)



    //N_pl phi/2 =N_particles+N_anyons/2
    //So phi =2*(N_particles+N_anyons/2)/N_pl
    size_t N_pl = (nw-1)*(nh-1);

    double phi =2.0*(N_particles+N_anyons*0.5)/N_pl;

    cout<<"Has phi "<<phi<<endl;

    //I strictly don't need to use uint64 here, but the compiler just gets real mad if I use different types
    for (uint64_t n = 0; n < N_states;++n)
    {
        //The hamiltonian is the sum of (t(j,k) a_j^+ a_k$ where the creation and annihilation operators give 0 unless a single one-distance swap can bring us from states n to state m, in which case it gives us t(j,k)
        state_t state = states[n];



        //Loop through all possible ways we can swap two particles <->
        for (int y = 0; y <  nh; ++y)
            for (int x = 0; x <  nw; ++x)
            {
                //Check for swaps <->
                if (x<nw-1)//No periodict boundary
                {
                    //destroy bit k and create bit j, here we only look at moving a bit <-, that only gives us the lower triangle n > m
                    int j = x+y*nw;
                    int k = x+1+y*nw;

                    //Binary math is uggly, but efficient
                    if (((state >> j) & ((state_t)1)) == 0 && ((state >> k) & ((state_t)1)) ==1)//By only checking this bit flip we stick to the lower triangle of the hamiltonian (this ensures that (original state number > new state number), the upper triangle come for free, as this is a hermitian matrix )
                    {
                        //do flip the bits, if this creates a new particle at j
                        //this is what the state would look like: uint64_t new_state = (state ^ (1UL << j)) ^ (1UL << k);
                        //Only look at lower triangle
                        //Now we want to find the index of this new state usign a binary search

                        int m=binSearch(states, (state ^ (  ((state_t)1) << j)) ^ (((state_t)1) << k), 0 ,N_states);
                                            //         ^^^^^^^^ ^      ^^^^^^^^^^
                                            //          |       |         |
                                            //Flip bit number j |         |
                                            //                  |         |
                                            //   Find bit number j   Flip bit number k
                        //Using t=1
                        complex number = -exp(-1i*M_PI*(2*phi)*double(y));//apparently c++ can not automatically cast int to double in this case

                        //We DESTROY at k and create at j, so the state with a bit at j is the left <- state, that state is state m, I believe the left state is the one which goes with the left index
                        H(m, n)+= number;
                        H(n,m)+= conj(number);
                    }
                }
                //Check for swaps
                // ^
                // |
                // v
                if (y<nh-1)
                {
                    //Swap bit i and j;
                    int j = x+y*nw;
                    int k = x+(y+1)*nw;

                    if (((state >> j) & ((state_t)1)) == 0 && ((state >> k) & ((state_t)1)) ==1)
                    {
                        //do flip the bits, if this creates a new particle at jActually, we don't need to save it
                        //uint64_t new_state = (state ^ (1UL << j)) ^ (1UL << k);
                        //Only look at lower triangle (original state number > new state number)
                        //Now we want to find the index of this new state usign a binary search
                        int m=binSearch(states, (state ^ ( ((state_t)1)<< j)) ^ ( ((state_t)1) << k), 0 ,N_states);//See previous swap for explanation

                        H(m, n)-= 1;
                        H(n,m)-= 1;//=conj(1)

                    }
                }


                //For the diagonal, check if we should add the potential of a particular site
                //Has a bit on this position
                int j = x+y*nw;
                if (((state >> j) & 1) == 1 )
                {
                    H(n, n)+= V(j);
                }
            }
        }


    //Somewhere to put the eigenvalues and vectors
    cx_vec eigval;
    cx_mat eigvec;

    //Here we get N_states eigenstates, we will only want to print the 3 lowest, so sort them in ascending order
    //eigs_gen = [eig]envalues of [s]parse [gen]eral (i.e. need not be real) matrix
    eigs_gen( eigval, eigvec, H, 1, "sr");//The lowest eigenstate


    this_eigval=eigval[0].real();
    this_eigvec = eigvec.col(0);

    //I don't know if the armadillo library guarantees that the eigenvectors are normalized, the documentation doesn't say so explicitly, as far as I can tell it always is, but I am not sure
    this_eigvec = normalise(this_eigvec );



    if ( !approx(eigval[0].imag(),0))
    {
        throw std::runtime_error("ERROR, energy has non-zero imaginary part");
    }

}
