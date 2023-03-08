using System;
using static System.Console;
using static System.Math;
using static matrix;
using static vector;

public static class minimizer
{
    //Mirror function to maximize instead
    public static (vector,int) max_qnewton(Func<vector,double>f, vector x0, double acc=1e-2, bool verbose=false,System.IO.StreamWriter writer =null)
    {
        Func<vector,double> F = (X)  => -f(X);
        return  qnewton(F,x0,acc,verbose,writer);
    }

    //Mirror function to maximize instead
    public static (vector,int) max_downhill_simplex(Func<vector,double>f, vector x0, double acc=1e-2, bool verbose=false,System.IO.StreamWriter writer =null,System.IO.StreamWriter writer_lowest =null)
    {
        Func<vector,double> F = (X)  => -f(X);
        return  downhill_simplex(F,x0,acc,verbose,writer,writer_lowest );
    }

    public static (vector,int) qnewton(Func<vector,double>f, vector x0, double acc=1e-2, bool verbose=false, System.IO.StreamWriter writer=null)
    {


        double little = Pow(2,-26);
        double alpha = 1e-4;//Used in acceptance condition for backtracking
        double epsilon = 1e-9;//Used to see if we should reset hessian

        vector x = x0.copy();//Better copy this, if the user wants to use the original input vector for something else.
        //Read the number of parameters the function outputs, and the number of parameters in the input

        int n = x.size;

        bool redo=true;//I have an emergency redo function, in case deltax is too small, in some rare cases, I have seen deltax be treated as a 0, breaking the algorithm uppon division, if that happen, we will retry here, with
        int max_steps=1000;//should never need that much

        int step = -1;//Start at -1 so we really start at 0
        while (redo)
        {
            redo=false;//In 99% of cases, this will only need to be done once, but maybe we are on a computer where machine epsilon is smaller than expected.


            step = -1;
            double fx=f(x);
            matrix B = new matrix(n,n);//Inverse hessian, start guess: identity

            //We will save the gradient between steps, as we need the gradient both at the start and at the end
            vector GradFx = new vector(n);

            if (writer != null)
                writer.WriteLine($"{x.asList()}\t{f(x)}");


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
                        vector offsetx = x.copy();
                        double deltax = Max(Abs(x[i])*little,0.25*little);//The tiny step to use when calculating the Jacobi matrix, note I do not enjoy 0 divisions, and it is possible that the user guesses that x[k]=0, so I add a floor to deltax

                        offsetx[i]+=deltax;

                        double foffsetxk=f(offsetx);
                        offsetx[i]+=deltax;
                        GradFx[i] = (foffsetxk-fx)/deltax;


                        if(double.IsInfinity(GradFx[i]))//This is very very rare
                        {
                            //CRAP, deltax was too small, and we got a 0 division, we must retry with a larger little number... unless we have done that too many times already

                            if (little>0.25)
                                WriteLine("ABORT: Got 0 division in Quasi Newton methods gradient too many times");

                            WriteLine("WARNING: Got 0 division in Quasi Newton method gradient, this may be caused by your computer having less precision for double point numbers than expected. Will retry with a larger Delta x");
                            redo = true;
                            x = resize(x0.copy(),n);//Reset our guess, and try again with a larger little number
                            little*=2;


                            break;
                        }
                    }

                if (!redo)//If finding the gradient worked
                {
                    if (verbose)
                    {
                        WriteLine($"Step {step}");
                        WriteLine(x.getString("x="));
                        WriteLine("");
                        WriteLine(GradFx.getString("Numeric ∇ f(x)="));
                        WriteLine("");
                    }

                    //Now need to do backtracking

                    // Let Delta x = - B ∇ f(x)
                    vector Dx = -B*GradFx;


                    // Let s = lambda Delta x
                    double lambda = 1;
                    vector S = Dx;
                    vector x_new = x+S;
                    bool linesearch_fail = true;

                    for (int i = 0; i < 32; ++i)
                    {
                        // Perform backtracking lambda -> lambda/2 and update s until f(x+s)<f(x)+alpha dot(s,∇ f(x)) (or until lambda too small)
                        if (f(x_new) < fx+ alpha*S.dot(GradFx)  )
                        {
                            linesearch_fail =false;
                            break;
                        }
                        lambda*=0.5;
                        S = lambda*Dx;
                        x_new = x+S;
                    }

                    x = x_new;
                    fx = f(x);
                    if (writer != null)
                        writer.WriteLine($"{x.asList()}\t{f(x)}");
                    if (verbose)
                    {
                        WriteLine("");
                        WriteLine(x.getString("This step x -> "));
                        WriteLine($"This step f(x) -> {fx}");
                    }

                    //With enough steps, this will actually work without resetting B. But we can do betterm now entering SR1 update


                    // CALCULATE ∇f(x + s) for the next step

                    vector GradFx_new = new vector(n);


                    for (int i = 0; i < n; ++i)
                    {
                        //Get the input, offset on the relevant element
                        vector offsetx = x.copy();
                        double deltax = Max(Abs(x[i])*little,0.25*little);//The tiny step to use when calculating the Jacobi matrix, note I do not enjoy 0 divisions, and it is possible that the user guesses that x[k]=0, so I add a floor to deltax

                        offsetx[i]+=deltax;

                        double foffsetxk=f(offsetx);
                        offsetx[i]+=deltax;
                        GradFx_new[i] = (foffsetxk-fx)/deltax;


                        if(double.IsInfinity(GradFx_new[i]))//This is very very rare
                        {
                            //CRAP, deltax was too small, and we got a 0 division, we must retry with a larger little number... unless we have done that too many times already

                            if (little>0.25)
                                WriteLine("ABORT: Got 0 division in Quasi Newton methods gradient too many times");

                            WriteLine("WARNING: Got 0 division in Quasi Newton method gradient, this may be caused by your computer having less precision for double point numbers than expected. Will retry with a larger Delta x");
                            redo = true;
                            x = resize(x0.copy(),n);//Reset our guess, and try again with a larger little number
                            little*=2;


                            break;
                        }
                    }


                    if (linesearch_fail)
                    {
                        B = new matrix(n,n);
                        if (verbose)
                        {
                            WriteLine("RESET B AS LINESEARCH FAILED");
                            WriteLine(B.getString("B="));
                        }

                    }
                    else
                    {

                    // Let y = ∇f(x + s)-∇f(x) and u = s- By
                    vector y = GradFx_new-GradFx;
                    vector u = S-B*y;


                    if (verbose)
                    {
                        WriteLine("");
                        WriteLine(GradFx_new.getString("∇f(x + s)"));
                        WriteLine("HAVE ");
                        WriteLine(y.getString("y="));
                        WriteLine("");
                        WriteLine(u.getString("u="));
                        WriteLine($"SO dot(y,u)={u.dot(y)}");
                    }

                    if (Abs(u.dot(y))>epsilon)
                    {

                        B = B+ matrix.outer_product(u,u)/u.dot(y);
                        if (verbose)
                        {
                            WriteLine("UPDATED B");
                            WriteLine(B.getString("B="));
                        }
                    }
                    else
                    {
                        if (verbose)
                        {
                            WriteLine("KEPT B");
                            WriteLine(B.getString("B="));
                        }
                    }

                    }

                    //Now, we already have whatever the gradient should be next step.
                    GradFx=GradFx_new;

                }


            if (GradFx.norm()<acc)
            {
                if (verbose)
                    WriteLine("BREAK DUE TO: gradient reached target");
                break;
            }

            if (step>=max_steps)//I am not sure how you would have > max steps ... ok actually you can, if we redo at the second to last step
            {
                if (verbose)
                    WriteLine("ALGORITHM FAILED: too many steps");
                break;
            }


            }//check if we should stop, from too many steps, or reaching the goal
            while (true );

        }

        if (writer != null)
            writer.Close();
        return (x,step);
    }

    public static (vector,int) downhill_simplex(Func<vector,double>f, vector x0, double acc=1e-2, bool verbose=false, System.IO.StreamWriter writer=null,System.IO.StreamWriter writer_lowest =null)
    {

        //Read the number of parameters the function outputs, and the number of parameters in the input
        int n = x0.size+1;
        int max_steps=1024;//should never need that much

        int step = -1;//Start at -1 so we really start at 0

        vector[] vertices = new vector[n+1];//The simplex vertices + the centroid
        double[] fs = new double[n];

        for (int i = 0; i<n ; ++i)
        {
            vertices[i] = x0.copy();


            if (i!=0)
            {
                vertices[i][i-1]+=1.0;//Start with an offset in each dimension, except for the first point
            }
            fs[i] = f(vertices[i]);

        }


        int lowest  = 0;//Which of the points are lowest and highest
        int highest = 0;

        //Check if any of the other actual points are higher/lower, and get the first centroid
        vertices[n] = vertices[0]/(n-1);//The centroid is (vertices[0]+...vertices[n-1])/(n-1) over all points which are NOT the highest
        for (int i = 1; i < n; ++i)
        {
            if (fs[i] < fs[lowest])
            {
                lowest = i;
            }
            if (fs[i] >= fs[highest])//The = here guarantees that the highest and lowest will not be the same, even if all points are identical
            {
                highest = i;
            }
            vertices[n] += vertices[i]/(n-1);//The centroid is (vertices[0]+...vertices[n-1])/(n-1) over all points which are NOT the highest
        }

        //Remove the highest point from the centroid sum
        vertices[n] -= vertices[highest]/(n-1);


        if (writer != null)
        {
            for (int i =0 ; i<n ; ++i)//Plot the triangle or whatever in a manner which gnuplot can understand
            {
                writer.WriteLine($"{vertices[i].asList()}\t{fs[i]}");
            }
            //For displaying this as a triangular grid with gnuplot
            //Return to start
            writer.WriteLine($"{vertices[0].asList()}\t{fs[0]}");
            //And then make sure we are not at the highest point
            writer.WriteLine($"{vertices[lowest].asList()}\t{fs[lowest]}");
        }
        if (writer_lowest != null)
            writer_lowest.WriteLine($"{vertices[lowest].asList()}\t{fs[lowest]}");




        bool do_end = false;
        bool reduced = false;//Was the last step reduction
        do
        {//Each step
            ++step;

            if (verbose)
            {
                WriteLine($"Step {step}:");
                for (int i =0 ; i<n ; ++i)//Plot the triangle or whatever in a manner which gnuplot can understand
                {
                    WriteLine($"{vertices[i].asList()}\t{fs[i]}"+(i == highest ?"\tHighest": ( (i==lowest) ? "\tLowest" : "")));
                }
            }

            //Try reflection
            vector new_vertex = vertices[n]-(vertices[highest]-vertices[n]);//This takes highest -> centroid - (vector from centroid to highest)
            double new_vertex_f = f(new_vertex);

            //Did it work
            if (new_vertex_f<fs[highest])
            {

                //If reflection works perfectly, try expansion
                if (new_vertex_f<fs[lowest])
                {//Oh we are going the right way, try expansion
                    lowest = highest;//Regardless of what happens, this is the lowest point now

                    //Try reflection
                    vector expanded= vertices[n]-2*(vertices[highest]-vertices[n]);//This takes highest -> centroid - 2*(vector from centroid to highest)
                    double expanded_f = f(expanded);

                    if (expanded_f<new_vertex_f)
                    {//accept expansion
                        new_vertex = expanded;
                        new_vertex_f = expanded_f;
                        if (verbose)
                            WriteLine($"Accept Expansion");

                    }
                    else
                        if (verbose)
                                WriteLine($"Reject Expansion, Accept Reflection");


                }
            }
            else//REJECT reflection
            {//Try contraction

                //Try contraction
                new_vertex = vertices[n]+0.5*(vertices[highest]-vertices[n]);//This takes highest -> centroid + 1/2 (vector from centroid to highest)
                new_vertex_f = f(new_vertex);

                if (new_vertex_f <fs[highest])
                {//accept expansion
                    if (verbose)
                        WriteLine($"Accept Contraction");

                    if (new_vertex_f <fs[lowest])
                    {
                        lowest = highest;//If, by change this new point is actually the lowest, reset the lowest marker to that
                    }

                }
                else
                {
                    if (verbose)
                        WriteLine($"Do reduction");

                    //Shrink everything towards the lowest point
                    for (int i = 0; i < n; ++i)
                    {
                        if (i != lowest)//The = here guarantees that the highest and lowest will not be the same, even if all points are identical
                        {
                            //This becomes the midpoint between this and the lowest
                            vertices[i] = 0.5*(vertices[i]+vertices[lowest]);



                        }
                    }

                    //also reset the centroid and the highest and lowest point, using a full recalculating as EVERYTHING changed
                    lowest  = 0;
                    highest = 0;

                    //Check if any of the other actual points are higher/lower, and get the first centroid
                    vertices[n] = vertices[0]/(n-1);//The centroid is (vertices[0]+...vertices[n-1])/(n-1) over all points which are NOT the highest

                    if ((vertices[0]-vertices[n-1]).norm()<acc)
                    {
                        if (verbose)
                            WriteLine($"Distance to neighbour within limit");
                        do_end=true;
                    }

                    for (int i = 1; i < n; ++i)
                    {
                        if (fs[i] < fs[lowest])
                        {
                            lowest = i;
                        }
                        if (fs[i] >= fs[highest])//The = here guarantees that the highest and lowest will not be the same, even if all points are identical
                        {
                            highest = i;
                        }
                        vertices[n] += vertices[i]/(n-1);//The centroid is (vertices[0]+...vertices[n-1])/(n-1) over all points which are NOT the highest

                        if ((vertices[i]-vertices[i-1]).norm()<acc)
                        {
                            if (verbose)
                                WriteLine($"Distance to neighbour within limit");
                            do_end=true;
                        }

                    }
                    //Remove the highest point from the centroid sum
                    vertices[n] -= vertices[highest]/(n-1);


                    reduced = true;



                }







            }

            if (!do_end && !reduced)//If we rand reduction, we already reset the lowest and highest, if we should end, we did not find a new vertex
            {//Reset the centroid and highest and lowest each step

                //If the new point is not the the lowest, the lowest is guaranteed to be unchanged

                vertices[highest] = new_vertex;
                //Reset highest and centroid
                fs[highest] =  new_vertex_f;

                //The highest may have changed, our first guess is that the new_vertex poitn is still the highest
                int new_vertex_point = highest;
                for (int i = 0; i < n; ++i)
                {
                    if (fs[i] >= fs[highest])//The = here guarantees that the highest and lowest will not be the same, even if all points are identical
                    {
                        highest = i;
                    }
                }

                if (new_vertex_point != highest)//Only change the centroid if the new point we added is not still the highest
                {
                    //Remove the point which now is the highest from the centroid sum
                    vertices[n] -= vertices[highest]/(n-1);
                    vertices[n] += vertices[new_vertex_point]/(n-1);
                }

                //Check if the point we just moved is too close to one neighbor (it is the only distance which changed, and it always has the same distance to both neighbour )

                int neighbour =(new_vertex_point == 0) ? n-1 : new_vertex_point -1;


                if ((vertices[new_vertex_point]-vertices[neighbour ]).norm()<acc)
                {
                    if (verbose)
                        WriteLine($"Distance to neighbour within limit");
                    do_end=true;
                }
            }




            if (writer != null)
            {
                for (int i =0 ; i<n ; ++i)//Plot the triangle or whatever in a manner which pyxplot can understand
                {
                    writer.WriteLine($"{vertices[i].asList()}\t{fs[i]}");
                }
                writer.WriteLine($"{vertices[0].asList()}\t{fs[0]}");//Close the triangle
                writer.WriteLine($"{vertices[lowest].asList()}\t{fs[lowest]}");//Then move to a point which will not be deleted next step, this will display a nice triangle grid
            }

            if (writer_lowest != null)
                writer_lowest.WriteLine($"{vertices[lowest].asList()}\t{fs[lowest]}");


            if (step>=max_steps)//I am not sure how you would have > max steps
            {
                if (verbose)
                    WriteLine("ALGORITHM FAILED: too many steps");
                do_end=true;
            }

        }
        while (!do_end);

        if (writer != null)
            writer.Close();
        return (vertices[lowest],step);
    }
}
