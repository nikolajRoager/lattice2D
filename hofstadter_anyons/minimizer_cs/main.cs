using System;
using static System.Console;
using static System.Math;
using static minimizer;
using static vector;

public static class main
{


    //double precision approximation
    public static bool approx(double a,double b,double tau=1e-9,double eps=1e-9)
    {
        if (Abs(a-b)<tau)
            return true;
        if (Abs(a-b)/(Abs(a)+Abs(b))<eps)
            return true;
        return false;
    }

    public static int Main(string[] args)
    {
        WriteLine($"----------------------------------------------");
        WriteLine($"Demonstrating Quasi newton minimization method");
        WriteLine($"----------------------------------------------");

        bool verbose = false;
        bool PASS = true;//Optimistic today are we?


        foreach(string s in args)
        {
            if (0==String.Compare(s,"-v"))
                verbose = true;
        }

        WriteLine($"Demonstration one, find max of f(x)=exp(-(x-2)^2) (true solution x=2");
        Func<vector,double> f0 = (X)  => Exp( -(X[0]-2)*(X[0]-2));
        WriteLine($"Starting at x0 = (0.0) with precision 10^-5. Now Running ...");
        (vector root0, int steps0) = max_qnewton(f0,new vector(0.0),1e-5,verbose );

        WriteLine("");
        WriteLine(root0.getString($"In {steps0} steps: Got predicted max at x="));

        if (root0.approx(new vector(2.0),1e-5,1e-5))//vector double is understood  as a vector with this element only
        {
            WriteLine("PASS this is within 10^-5 of 2.0");
        }
        else
        {
            WriteLine("FAIL this is not within 10^-5 of 2.0");
            PASS=false;
        }

        WriteLine($"----------------------------------------------");

        WriteLine($"Demonstration Rosenbrock, find minimum of f(x,y)=((1-x)^2+100*(y-x^2)^2) (minimum is at 1,1)");
        Func<vector,double> f1 = (X)  => ((1-X[0])*(1-X[0])+100*(X[1]-X[0]*X[0])*(X[1]-X[0]*X[0]) ) ;
        WriteLine($"Starting at x0 = (-1,2) with precision 10^-5. Now Running ...");
        (vector root1, int steps1) = qnewton(f1,new vector(-1,2),1e-5,verbose, new System.IO.StreamWriter("rosenbrock.tsv"));

        WriteLine("");
        WriteLine(root1.getString("Got predicted root x="));

        WriteLine("");
        WriteLine($"In {steps1} steps: Has f(x)={f1(root1)}");
        WriteLine("");
        if (root1.approx(new vector(1.0,1.0),1e-5,1e-5))
        {
            WriteLine("PASS this is within 10^-5 of (1,1)");
        }
        else
        {
            WriteLine("FAIL this is not within 10^-5 of (1,1)");
            PASS=false;
        }



        WriteLine($"----------------------------------------------");
        WriteLine($"Demonstration Himmelblau, find minimum of f(x,y)=(x^2+y-11)^2+(x+y^2-7)^2\nminima at: (3.0,2.0),(-2.805118, 3.131312),(-3.779310, -3.283186) and (3.584428, -1.848126)");
        Func<vector,double> f2 = (X)  => ((X[0]*X[0]+X[1]-11)*(X[0]*X[0]+X[1]-11)+(X[0]+X[1]*X[1]-7)*(X[0]+X[1]*X[1]-7));
        WriteLine($"Starting at x0 = (0,0) with precision 10^-5. Now Running ...");
        (vector root2, int steps2) = qnewton(f2,new vector(0,0),1e-5,verbose , new System.IO.StreamWriter("himmelblau.tsv"));

        WriteLine("");
        WriteLine(root2.getString("Got predicted root x="));

        WriteLine("");
        WriteLine($"In {steps2} steps: Has f(x)={f2(root2)}");
        WriteLine("");
        if ( approx(f2(root2),0.0,1e-5,1e-5))
        {
            WriteLine("PASS the function is within 10^-5 of 0, here; which is the known minimum");
        }
        else
        {
            PASS=false;
            WriteLine("FAIL the function is not within 10^-5 of 0, here; which is the known minimum");
        }
        WriteLine($"----------------------------------------------");
        if (PASS)
            WriteLine("ALL TESTS PASSED");
        else
            WriteLine("SOME TESTS FAILED");

        return 0;
    }

}
