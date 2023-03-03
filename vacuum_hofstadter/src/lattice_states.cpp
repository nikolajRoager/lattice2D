#include"generate_states.hpp"

//Has already been done in my header, but I like to include everything explicitly in every file (the compiler will ignore it)
#include<cstdint>

//NOTE! std::vector is a dynamic sized list, it is not a vector in the linear algebra sense of the word, that will be the class vec from the armadillo library
#include<vector>
#include<iostream>
#include<string>
#include<complex>

//The actual math library
#include <armadillo>

using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)


bool approx(double a, double b,double tau = 1e-9, double epsilon =1e-9)
{

    if (abs(a-b)<=tau)
        return true;
    else if (abs(a-b)/(abs(a)+abs(b))<=epsilon)//No risk of divide by 0, if the user asks for a=b=0, then the above always returns true (unless some idiot asks for a negative tolerance and then asks if 0=0)
        return true;
    else
        return false;
}

int main(int argc, char* argv[])
{
    //Read arguments
    if (argc!=4 && argc!=5)
    {
        cout<<"Invalid, argument, need "<<argv[0]<<" width height M (Verbose)"<<endl;
        return 1;
    }
    bool verbose = false;
    if (argc==5)
    {
        if (string ("verbose").compare(argv[4]) == 0)
            verbose=true;
    }

    //If this fails, we get 0
    int nw = atoi(argv[1]);
    int nh = atoi(argv[2]);

    int N_sites = nh*nw;
    int N_particles = atoi(argv[3]);


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

    //I am very unlikely to get a number of states which needs uint64_t here, but I do need to use the same type as I use for the states (written as a 64 bit number), some of the binary operations may fail otherwise
    uint64_t N_states=0;//Will be generated
    vector<state_t> states(N_states);

    generate(states,N_states,N_sites,N_particles);

    cout <<"Found "<<N_states<<" states"<<endl;



    if (verbose)
    {
        cout<<"\nAll states:"<<endl;
        //I strictly don't need to use uint64 here, but the compiler just gets real mad if I use different types
        for (uint64_t i = 0; i < N_states; ++i)
        {
            print_state(states[i],nw,nh);
        }

    }

    cout<<"\n\n Now creating hermitian hamiltonian ("<<N_states<<" by "<<" "<<N_states<<")"<<endl;

    //Make our hamiltonian matrix
    sp_cx_mat H (N_states, N_states);//Sparse matrix, as we have mostly 0 entries, initialized as only zeros by default

    //Potential
    vec V(N_sites);
    //Loop through all our states and calculate all matrix elements
    //In my case, I know only nearest neighbours interact (Hofstadter model), so I only loop through each state once, and I only consider neighbour with a lower state than this state (Lower triangle, we know the matrix is hermitian! so we will just mirror everything)



    //N_pl phi/2 =N_particles+S/2
    //So phi =2*(N_particles+S/2)/N_pl
    size_t N_pl = (nw-1)*(nh-1);

    size_t S=0;

    double phi =2.0*(N_particles+S*0.5)/N_pl;

    cout<<"Has φ="<<phi<<endl;


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

    // (compare to manual calculations), print the hamiltonian
    if (verbose)
    {
        cout<<"Hamiltonian:"<<endl;
        cout<<H<<endl;
    }


    printf("\n Now diagonalizing hamiltonian...\n");


    //Somewhere to put the eigenvalues and vectors
    cx_vec eigval;
    cx_mat eigvec;

    //Here we get N_states eigenstates, we will only want to print the 3 lowest, so sort them in ascending order
    //eigs_gen = [eig]envalues of [s]parse [gen]eral (i.e. need not be real) matrix
    eigs_gen( eigval, eigvec, H, 3, "sr");//3,sr= 3 eigenvalues with smallest real part of the eigenvalues, hopefully the eigenvalues are REAL, but the function does not know that




    cout << "\n ... Done, here are the 3 lowest energies and eigenvectors"<<endl;


    for (int i = 0; i < 3; ++i)
    {
        cout<<"E_"<<i<<"="<<eigval[i].real()<<endl;


        cx_vec this_eigvec = eigvec.col(i);

        //I don't know if the armadillo library guarantees that the eigenvectors are normalized, the documentation doesn't sasy so, as far as I can tell it always is, but I am not sure
        this_eigvec = normalise(this_eigvec );


        vector<double> prop_list(N_states);



        //Verify that this is an eigenvector
        cx_vec shouldBe0 = H*this_eigvec-eigval[i].real()*this_eigvec;//Use the real part of the energy, because the imaginary part is naught but a numerical error


        if ( !approx(eigval[i].imag(),0))
        {
            cout<<"ERROR, energy has non-zero imaginary part"<<endl;
            return 1;
        }


        if (verbose)
        {
            cout<<"Corresponding Eigenstate "<<endl;
            cout<<this_eigvec<<endl;
        }


        double density_sum=0;
        //I strictly don't need to use uint64 here, but the compiler just gets real mad if I use different types
        for (uint64_t j = 0; j < N_states; ++j)
        {
            prop_list[j]= norm(this_eigvec[j]);

            if ( !approx(shouldBe0[j].real(),0)|| !approx(shouldBe0[j].imag(),0))
            {
                cout<<"ERROR, the eigenstate found is NOT an eigenstate"<<endl;
                return 1;
            }

            density_sum+=prop_list[j]*N_particles;
        }

        cout<<"Total particle density is "<<density_sum<<endl;
        cout<<"All test PASSED, is valid eigenstate and eigenvalue"<<endl;


        //Save the first three states to stderr (it is not actually an error, but this is a way I can print this data to another file, which I can plot)
        print_superpositions_data(states,prop_list , N_states,  nw, nh);
        cerr << endl;
        cout << endl;
    }

    return 0;
}
