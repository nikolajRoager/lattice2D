#include "minimize.hpp"

#include<cstdint>
#include<iostream>
#include<string>
#include<armadillo>
#include<functional>

using state = uint64_t;
using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)


vec qnewton( function<double(vec)> func , vec x0, size_t& step, double acc, bool verbose)
{
    double little = pow(2,-26);
    double alpha = 1e-4;//Used in acceptance condition for backtracking
    double epsilon = 1e-9;//Used to see if we should reset hessian

    vec x = x0;//Calls copy constructor by default

    //Read the number of parameters the function outputs, and the number of parameters in the input
    int n = x.size();


    bool redo=true;//I have an emergency redo function, in case deltax is too small, in some rare cases, I have seen deltax be treated as a 0, breaking the algorithm uppon division, if that happen, we will retry here
    int max_steps=256;//should never need that much

    while (redo)
    {
        redo=false;//In 99% of cases, this will only need to be done once, but maybe we are on a computer where machine epsilon is smaller than expected.


        step = -1;
        double fx=func(x);
        mat B = mat(n,n,fill::eye);//Inverse hessian, start guess: identity

        //We will save the gradient between steps, as we need the gradient both at the start and at the end
        vec GradFx = vec(n);


        if (verbose)
        {
            //cout<<" x =";
            //for (size_t i = 0; i<n;++i)
            //    cout<<' '<<x[i];

            cout<<"Step "<<step<<"\t f(x)="<<func(x)<<endl;
        }
        do
        {//Each step
            ++step;
            //PART 1: Calculate the Gradient
            // ∇ f(x)
            //
            // PART 2: backtracking:
            // Let Delta x = - B ∇ f(x)
            // Let s = lambda Delta x
            // Perform backtracking lambda -> lambda/2 and update s until f(x+s)<f(x)+alpha dot(s,∇ f(x)) (or until lambda too small)
            //
            // Symmetric Rank 1 update:
            //
            // CALCULATE ∇f(x + s) again !!?
            //
            // Let y = ∇f(x + s)-∇f(x) and u = s- By
            //
            // if dot(u,y)>epsilon
            //     let B = outer_product(u,u)/dot(u,y)
            // else
            //     let B = 1
            // endif


            //If this is the first step, we need to get gradient here matrix
            if  (step==0)
                for (int i = 0; i < n; ++i)
                {
                    //Get the input, offset on the relevant element
                    vec offsetx = x;
                    double deltax = std::max(std::abs(x[i])*little,0.25*little);//The tiny step to use when calculating the Jacobi matrix, note I do not enjoy 0 divisions, and it is possible that the user guesses that x[k]=0, so I add a floor to deltax

                    offsetx[i]+=deltax;

                    double foffsetxk=func(offsetx);
                    offsetx[i]+=deltax;
                    GradFx[i] = (foffsetxk-fx)/deltax;

                    if(std::isinf(GradFx[i]))//This is very very rare
                    {
                        //CRAP, deltax was too small, and we got a 0 division, we must retry with a larger little number... unless we have done that too many times already

                        if (little>0.25)
                            cout<<("ABORT: Got 0 division in Quasi Newton methods gradient too many times")<<endl;

                        cout<<("WARNING: Got 0 division in Quasi Newton method gradient, this may be caused by your computer having less precision for double point numbers than expected. Will retry with a larger Delta x")<<endl;
                        redo = true;
                        x = x0;//Reset our guess, and try again with a larger little number
                        little*=2;


                        break;
                    }
                }

            if (!redo)//If finding the gradient worked
            {
                if (verbose)
                {
                    //cout<<"Step "<<step<<endl;
                    //cout<<"x= "<<x<<endl<<endl;
                    //cout<<"numeric = ∇ f(x)="<<GradFx<<endl<<endl<<endl;
                }

                //Now need to do backtracking

                // Let Delta x = - B ∇ f(x)
                vec Dx = -B*GradFx;


                // Let s = lambda Delta x
                double lambda = 1;
                vec S = Dx;
                vec x_new = x+S;


                bool linesearch_fail = true;

                for (int i = 0; i < 32; ++i)
                {
                    // Perform backtracking lambda -> lambda/2 and update s until f(x+s)<f(x)+alpha dot(s,∇ f(x)) (or until lambda too small)
                    if (func(x_new) < fx+ alpha*dot(S,GradFx)  )
                    {
                        linesearch_fail =false;
                        break;
                    }
                    lambda*=0.5;
                    S = lambda*Dx;
                    x_new = x+S;
                }

                x = x_new;
                fx = func(x);
                //if (writer != null)
                //    writer.WriteLine($"{x.asList()}\t{f(x)}");
                if (verbose)
                {
                    cout<<"Step "<<step<<"\t f(x)="<<func(x)<<endl;
                    //cout<<endl;
                    //cout <<"This step x -> "<<x<<endl;
                    //cout<<"This step f(x) -> "<<fx<<endl;
                }

                //With enough steps, this will actually work without resetting B. But we can do betterm now entering SR1 update


                // CALCULATE ∇f(x + s) for the next step

                vec GradFx_new = vec(n);


                for (int i = 0; i < n; ++i)
                {
                    //Get the input, offset on the relevant element
                    vec offsetx = x;
                    double deltax = std::max(std::abs(x[i])*little,0.25*little);//The tiny step to use when calculating the Jacobi matrix, note I do not enjoy 0 divisions, and it is possible that the user guesses that x[k]=0, so I add a floor to deltax

                    offsetx[i]+=deltax;

                    double foffsetxk=func(offsetx);
                    offsetx[i]+=deltax;
                    GradFx_new[i] = (foffsetxk-fx)/deltax;


                    if(std::isinf(GradFx_new[i]))//This is very very rare
                    {
                        //CRAP, deltax was too small, and we got a 0 division, we must retry with a larger little number... unless we have done that too many times already

                        if (little>0.25)
                            cout<<"ABORT: Got 0 division in Quasi Newton methods gradient too many times"<<endl;

                        cout<<"WARNING: Got 0 division in Quasi Newton method gradient, this may be caused by your computer having less precision for double point numbers than expected. Will retry with a larger Delta x"<<endl;
                        redo = true;
                        x = x0;//Reset our guess, and try again with a larger little number
                        little*=2;


                        break;
                    }
                }


                if (linesearch_fail)
                {
                    B = mat(n,n,fill::eye);
                    if (verbose)
                    {
                        cout<<"RESET B AS LINESEARCH FAILED"<<endl;
                        //cout<<"B="<<B<<endl;
                    }

                }
                else
                {

                    // Let y = ∇f(x + s)-∇f(x) and u = s- By
                    vec y = GradFx_new-GradFx;
                    vec u = S-B*y;


                    if (verbose)
                    {
                        //cout<<endl;
                        //cout<<"∇f(x + s)="<< GradFx_new<<endl;
                        //cout<<"y="<<y<<endl;
                        //cout<<endl;
                        //cout<<"u="<<u<<endl;
                        //cout<<"So dot(y,u)={u.dot(y)}"<<endl;
                    }

                    if (std::abs(dot(u,y))>epsilon)
                    {

                        B = B+ mat(u)*trans(mat(u))/dot(u,y);
                        if (verbose)
                        {
                            //cout<<"UPDATED B"<<endl;
                            //cout<<B<<endl;
                        }
                    }
                    else
                    {
                        if (verbose)
                        {
                            //cout<<"KEPT B"<<endl;
                            //cout<<B<<endl;
                        }
                    }

                }

                //Now, we already have whatever the gradient should be next step.
                GradFx=GradFx_new;

            }


        if (norm(GradFx)<acc)
        {
            if (verbose)
                cout<<"BREAK DUE TO: gradient reached target"<<endl;
            break;
        }

        if (step>=max_steps)//I am not sure how you would have > max steps ... ok actually you can, if we redo at the second to last step
        {
            if (verbose)
                cout<<"ALGORITHM FAILED: too many steps"<<endl;
            break;
        }

        }//check if we should stop, from too many steps, or reaching the goal
        while (true);

    }

    return x;
}



//Mirror function to maximize instead
vec max_qnewton( function<double(vec)> func , vec x0, size_t& steps, double acc, bool verbose)
{
    return  qnewton( [func](vec x) -> double {return -func(x);},x0,steps,acc,verbose);
}

//Versions without counting number of steps
vec qnewton( function<double(vec)> func , vec x0, double acc, bool verbose)
{
    size_t steps;
    return  qnewton( func,x0,steps,acc,verbose);
}

//Mirror function to maximize instead
vec max_qnewton( function<double(vec)> func , vec x0, double acc, bool verbose)
{
    size_t steps;
    return  qnewton( [func](vec x) -> double {return -func(x);},x0,steps,acc,verbose);
}
