#include<generate_states.hpp>

//Has already been done in my header, but I like to include everything explicitly in every file (the compiler will ignore it)
#include<cstdint>
#include<vector>
#include<iostream>

using state = uint64_t;
using namespace std;

int main(int argc, char* argv[])
{
    //Read arguments
    if (argc!=4)
    {
        cout<<"Invalid, argument, need "<<argv[0]<<" width height N"<<endl;
        return 1;
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
    if (N_sites >64)
    {
        cout <<"This program only allows up to 64 sites, got  "<<nw<<"x"<<nh<<"="<<N_sites<<endl;
        return 1;
    }

    if (N_sites < N_particles)
    {
        printf("More particles %i than lattice sites  %ix%i=%i\n",N_particles,nw,nh,N_sites);
        return 1;
    }


    if (0 > N_particles)
    {
        printf("Negative particle number %i not allowed",N_particles);
        return 1;
    }

    cout << nw<<"x"<<nh<<" Lattice with "<<N_particles<<" particles for "<<N_sites<<" sites"<<endl;

    cout << "Generating and printing ALL legal states to stderr"<<endl;

    uint64_t N_states=0;//Will be generated
    vector<state> states(N_states);

    generate(states,N_states,N_sites,N_particles);

    cout <<"Found "<<N_states<<" states"<<endl;


    for (int i = 0; i < N_states; ++i)
    {
        print_state_data(states[i],nw,nh);
        cerr<<endl;

        //Count that we always have the correct number of particles
        int particles=0;
        for ( uint64_t l = 1lu; l!=0; l=l<<1)
            if (states[i] & l)
                ++particles;

        if (particles!=N_particles)
        {
            cout<<"ERROR Wrong number of particles, expected "<<N_particles<<" counted "<<particles<<endl;
            return 1;
        }



    }

    return 0;
}
