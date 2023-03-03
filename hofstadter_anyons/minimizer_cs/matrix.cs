using System;
using static System.Math;
//using static System.Console;

public class matrix
{
    public readonly int height,width;
    private double[] data;

    //Generate rectangular matrix
    public matrix(int n, int m)
    {
		height=n;
        width=m;
		data = new double[height*width];

        //Default to identity
        for (int i = 0; i < height; ++i)
            for (int j = 0; j < width; ++j)
                this[i,j]= (i==j) ? 1 : 0;
    }


    public matrix(int n)
    {
		height=n;
        width=n;
		data = new double[height*width];

        //Default to identity
        for (int i = 0; i < height; ++i)
            for (int j = 0; j < width; ++j)
                this[i,j]= (i==j) ? 1 : 0;
    }

    //Create get and set functions
	public double this[int i,int j]
    {
		get => data[i+j*height];
		set => data[i+j*height]=value;
	}

    //Just flipping the arguments will also do
    public matrix transpose()
    {
        matrix Out = new matrix(width,height);

        for (int i = 0; i < height; ++i)
            for (int j = 0; j < width; ++j)
                Out[j,i]=this[i,j];
        return Out;
    }
    //Pretty display function
    public override string ToString()
    {
        string Out = "";
        for (int i = 0; i < height; ++i)
        {
            Out+= "  |";
            for (int j = 0; j < width; ++j)
            {
                Out+=((Abs(this[i,j])< 100) ? " " : "")+((Abs(this[i,j])< 10) ? " " : "")+((this[i,j]>=0) ? " " : "")+string.Format(" {0:N5}",this[i,j]);
            }
            if (i<height-1)
                Out+=" |\n";
            else
                Out+= " |";
        }
        return Out;
    }

    //Pretty display function, with something in front, i.e. A = ...
    public string getString(string prefix)
    {
        string Out = "";
        for (int i = 0; i < height; ++i)
        {

            if (i != height/2)
                for (int j = 0; j<prefix.Length; ++j)
                    Out+= " ";
            else
                Out+=prefix;

            Out+= "|";
            for (int j = 0; j < width; ++j)
            {

                Out+=((Abs(this[i,j])< 100) ? " " : "")+((Abs(this[i,j])< 10) ? " " : "")+((this[i,j]>=-0.001) ? " " : "")+string.Format(" {0:N5}",this[i,j]);
            }
            if (i<height-1)
                Out+=" |\n";
            else
                Out+= " |";
        }
        return Out;
    }


    //Same as above, but returns a list of strings in order to display matrix multiplication pretty
    public string[] getStrings(string prefix)
    {
        string[] Out = new string[height];
        for (int i = 0; i < height; ++i)
        {
            Out[i]="";
            if (i != height/2)
                for (int j = 0; j<prefix.Length; ++j)
                    Out[i]+= " ";
            else
                Out[i]+=prefix;

            Out[i]+= "|";
            for (int j = 0; j < width; ++j)
            {

                Out[i]+=((Abs(this[i,j])< 100) ? " " : "")+((Abs(this[i,j])< 10) ? " " : "")+((this[i,j]>=-0.001) ? " " : "")+string.Format(" {0:N5}",this[i,j]);
            }
            Out[i]+=" |";
        }
        return Out;
    }


    //Python compatible display version, for testing
    public string ToPython()
    {
        //Python compatible display version, for testing
        string Out = "np.matrix([";
        for (int i = 0; i < height; ++i)
        {
            Out+= "[";
            for (int j = 0; j < width; ++j)
            {

                Out+=string.Format(" {0}",this[i,j]);
                if (j<width-1)
                    Out+=",";
            }
            if (i<height-1)
                Out+=" ],";
            else
                Out+= " ]])";
        }
        return Out;
    }

    //Pretty display of the linear equation
    public string getString(matrix b)
    {
        string Out = "";
        if (b.width!=1)
            throw new ArgumentException("b should be a single column");
        if (b.height!=width)
            throw new ArgumentException("b should have same height as matrix width");
        for (int i = 0; i < height; ++i)
        {


            Out+= "   ";
            for (int j = 0; j < width; ++j)
            {

                Out+=((this[i,j]>=0 && j>0) ? "+" : "")+string.Format(" {0:N5}",this[i,j])+$" * x[{j}] ";
            }
            if (i<height-1)
                Out+=$" = {b[i,0]} \n";
            else
                Out+= $" = {b[i,0]}";
        }
        return Out;
    }


    public void randomize_symmetric()
    {
        var generator = new Random();

        if (width != height)
            throw new ArgumentException("Symmetric matrices must be square");

        for (int j = 0; j < width; ++j)
            for (int i = j; i < width; ++i)
            {
                this[i,j]=generator.NextDouble();
                this[j,i]=this[i,j];
            }
    }

    public void randomize()
    {
        var generator = new Random();

        for (int j = 0; j < width; ++j)
            for (int i = 0; i < height; ++i)
                this[i,j]=generator.NextDouble();
    }

    public double[] get_data ()
    {
        return data;
    }


    public double colNorm (int j)
    {
        double Norm2 =0;


        for (int i = 0; i < height; ++i)
        {
            Norm2+= this[i,j]*this[i,j];

        }
        return Sqrt(Norm2);
    }

    //Essentially the copy constructor which would be used by default in C++
    public matrix copy()
    {
        matrix Out = new matrix(height,width);
        for(int i=0;i<height;i++)
            for(int j=0;j<width;j++)
                Out[i,j]=this[i,j];
        return Out;
    }

    //Run Gram Schmidt on this matrix
    public (matrix,matrix) getQR()
    {

        //The algorithm assumes that the columns of A are linearly independent, lets just hope they are ...

        matrix Q= this.copy();
        matrix R= new matrix(width,width);

        if (height<width)
            throw new ArgumentException("matrix height must be greater or equal to width for this implementation to work");


        for (int i = 0; i < width; ++i)
        {
            R[i,i]=Q.colNorm(i);
            //Console.WriteLine($"R[{i},{i}] :{R[i,i]}");

            //Set column vector i in Q to be A normalized
            for (int k = 0; k < height; ++k)
            {
                Q[k,i]=Q[k,i]/R[i,i];
            }
            for (int j = i+1; j < width; ++j)
            {

                R[i,j]=0;
                for (int k = 0; k < height; ++k)
                {
                    R[i,j]+=Q[k,i]*Q[k,j];
                }
                for (int k = 0; k < height; ++k)
                {
                    Q[k,j]-=Q[k,i]*R[i,j];
                }
            }
        }

        return (Q,R);
    }

    public static vector QRsolve (matrix Q, matrix R, vector b)
    {
        if (b.size!=Q.height)
            throw new ArgumentException("b should have same height as Q matrix");
        if (Q.width!=R.width)
            throw new ArgumentException("Q and R should have same width");


        //Already assumes that Q is orthagonal, if not we are just solving   R*x = Q^T b
        var c = Q.transpose()*b;
        var x = c.copy();
        //Now R x = c can be solved with back-substitution
        //As in the chapter, I use in-place back substitution

        //Error.WriteLine("\n");
        //Error.WriteLine(c.getString("DO TRY TO SOLVE R*x = Q^T b ="));
        for (int i = c.size-1; i>=0; --i)
        {
            double sum = 0;
            for (int k = i+1; k<c.size; ++k)
                sum+=R[i,k]*x[k];
            x[i]=(x[i]-sum)/R[i,i];
        }

        return x;
	}


    public static matrix inverse(matrix Q, matrix R)
    {
        if (Q.width!=R.width)
            throw new ArgumentException("Q and R should have same width");

        //We really want to solve $height equations QR x = [0,0,...1 ... 0]

        //We can write this as a matrix equation, which can be solved with back substitution
        matrix unity = new matrix(Q.height,Q.width);
        //Already assumes that Q is orthagonal
        //Now R x = c can be solved with back-substitution on the entire matrix
        matrix Out = Q.transpose()*unity;
        for (int j = 0; j < Q.width; ++j)
        {
            for (int i = Out.height-1; i>=0; --i)
            {
                double sum = 0;
                for (int k = i+1; k<Out.height; ++k)
                    sum+=R[i,k]*Out[k,j];
                Out[i,j]=(Out[i,j]-sum)/R[i,i];
            }
        }
        return Out;
	}

    //I did not note down the definition pseudo-inverse's, but it sounds somewhat similar to penrose inverse
    public static matrix penrose_inverse(matrix Q, matrix R)
    {
        //We have do R^T R A^p = (A)^T

        //R^T R A^p = (Q R)^T
        matrix RT = R.transpose();
        matrix AT = (Q*R).transpose();



        //We can write this as a matrix equation, which can be solved with forward+back substitution

        matrix OUT = AT.copy();

        //Forwards substitution,column by column, in place
        for (int j = 0; j < OUT.width; ++j)
        {
            for (int i = 0; i< OUT.height; ++i)
            {
                double sum = 0;
                for (int k = 0; k<i; ++k)
                    sum+=RT[i,k]*OUT[k,j];
                OUT[i,j]=(OUT[i,j]-sum)/RT[i,i];
            }
        }
/*
        matrix TEST = RT*OUT;
        Error.WriteLine("After forward substitution now OUT = R A^p");
        Error.WriteLine(TEST.getString("RT*OUT   ="));
        Error.WriteLine("\n");
        Error.WriteLine(AT.getString("A^T (same) ="));
*/

        //Now we have: R A^p = (R A^p), use backwards substitution to get A^P
        //Backward substitution in place
        for (int j = 0; j < OUT.width; ++j)
        {
            for (int i = OUT.height-1; i>=0; --i)
            {
                double sum = 0;
                for (int k = i+1; k<OUT.height; ++k)
                    sum+=R[i,k]*OUT[k,j];
                OUT[i,j]=(OUT[i,j]-sum)/R[i,i];
            }
        }
    /*
        TEST = RT*R*OUT;
        Error.WriteLine("After backwards substitution now OUT = A^p");
        Error.WriteLine(TEST.getString("R^T*R*OUT ="));
        Error.WriteLine("\n");

        Error.WriteLine((AT).getString("A^T (same)="));
        Error.WriteLine("\n");
*/

        return OUT;
	}

    //Not as detailed operations as in the matlib class, but enough to get by:
    public static matrix operator+(matrix a, matrix b)
    {
        matrix Out = new matrix(a.height,a.width);
        if (a.height!=b.height || a.width!=b.width)
            throw new ArgumentException("matrix height width do not match for addition");
        for(int i=0;i<a.height;i++)
            for(int j=0;j<a.width;j++)
                Out[i,j]=a[i,j]+b[i,j];
        return Out;
	}

    public static matrix operator-(matrix a)
    {
        matrix Out = new matrix(a.height,a.width);

        for(int i=0;i<a.height;i++)
            for(int j=0;j<a.width;j++)
                Out[i,j]=-a[i,j];
        return Out;
	}

    public static matrix operator- (matrix a, matrix b)
    {
        matrix Out = new matrix(a.height,a.width);
        if (a.height!=b.height || a.width!=b.width)
            throw new ArgumentException("matrix height width do not match for subtraction");

        for(int i=0;i<a.height;i++)
            for(int j=0;j<a.width;j++)
                Out[i,j]=a[i,j]-b[i,j];
        return Out;
	}

    //Matrix scaling
    public static matrix operator*(matrix a, double x)
    {
        matrix Out = new matrix(a.height,a.width);
        for(int i=0;i<a.height;i++)
            for(int j=0;j<a.width;j++)
                Out[i,j]=a[i,j]*x;
        return Out;
    }
    //Matrix scaling
    public static matrix operator/(matrix a, double x)
    {
        matrix Out = new matrix(a.height,a.width);
        for(int i=0;i<a.height;i++)
            for(int j=0;j<a.width;j++)
                Out[i,j]=a[i,j]/x;
        return Out;
    }

    //Is this matrix evil (NaN or infinity)
    public bool isEvil()
    {
        for(int i=0;i<height;i++)
            for(int j=0;j<width;j++)
                if (double.IsInfinity(this[i,j]) || double.IsNaN(this[i,j]))
                    return true;
        return false;
    }

    //same idea as in the matrix library given, allow multiplication from both sides
    public static matrix operator*(double x, matrix a){ return a*x; }

    //Matrix Matrix multiplication, this is the important part
    public static matrix operator* (matrix a, matrix b)
    {
        //Need mathcing matrices
        if (a.width!=b.height )
            throw new ArgumentException($"matrices need to have matching width/height for multiplication to be defined, got {a.height} by {a.width} times {b.height} by {b.width}");


        matrix Out = new matrix(a.height,b.width);
        for (int i=0;i<a.height;i++)
            for (int j=0;j<b.width;j++)
            {
                Out[i,j]=0;
                for (int k=0;k<a.width;k++)
                    Out[i,j]+=a[i,k]*b[k,j];
            }

        return Out;
    }

    //Outer vector product, creates a matrix
    public static matrix outer_product (vector a, vector b)
    {
        //Here we assume a is a column and b is a row
        //Need mathcing matrices

        matrix Out = new matrix(a.size,b.size);
        for (int i=0;i<a.size;i++)
            for (int j=0;j<b.size;j++)
            {
                //a.width = 1, so only this
                Out[i,j]=a[i]*b[j];
            }

        return Out;
    }



    //Multiply matrix on (column) vectors
    public static vector operator* (matrix a, vector b)
    {
        //Need mathcing matrices
        if (a.width!=b.size )
            throw new ArgumentException($"Matrix needs same width as column vector size got {a.height} by {a.width} times {b.size}");


        vector Out = new vector(a.height);
        for (int i=0;i<a.height;i++)
            {
                Out[i]=0;
                for (int k=0;k<a.width;k++)
                    Out[i]+=a[i,k]*b[k];
            }

        return Out;
    }

    //Multiply (row) vectors on a matrix
    public static vector operator* (vector a, matrix b)
    {
        //Need mathcing matrices
        if (a.size!=b.height)
            throw new ArgumentException($"Matrix needs same height as row vector got {b.height} by {b.width} times {a.size}");

        vector Out = new vector(b.width);
        for (int j=0;j<b.width;j++)
        {
            Out[j]=0;
            for (int k=0;k<a.size;k++)
                Out[j]+=a[k]*b[k,j];
        }


        return Out;
    }


    //double precision approximation
    public static bool approx(double a,double b,double tau=1e-5,double eps=1e-5)
    {
        if (Abs(a-b)<tau)
            return true;
        if (Abs(a-b)/(Abs(a)+Abs(b))<eps)
            return true;
        return false;
    }

    public bool approx(matrix that)
    {
        if (height != that.height || width != that.width)
            return false;
        for (int i = 0; i < height; ++i)
            for (int j = 0; j < width; ++j)
                if (!approx(that[i,j],this[i,j]))
                    return false;
        return true;
    }


    public bool approx(matrix A,matrix B)
    {
        if (A.height != B.height || B.width != A.width)
            return false;
        for (int i = 0; i < A.height; ++i)
            for (int j = 0; j < A.width; ++j)
                if (!approx(A[i,j],B[i,j]))
                    return false;
        return true;
    }


    public bool is_uptriangle()
    {
        for (int i = 0; i < height; ++i)
            for (int j = 0; j < i; ++j)
                if (!approx(this[i,j],0))
                    return false;

        return true;
    }


}
