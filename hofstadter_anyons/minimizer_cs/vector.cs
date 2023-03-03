using System;
using static System.Math;

//A dumped down version of the matrix class with only one column, you could just use a column matrix for all of this
public class vector
{
    public readonly int size;
    private double[] data;




    //Read vector from an input array
    public vector(params double[] _data)
    {
		size=_data.Length;

		data = new double[size];
        for (int i = 0; i < size; ++i)
            this[i]= _data[i];
    }
    public vector(int n)
    {
		size=n;
		data = new double[size];

        for (int i = 0; i < size; ++i)
            this[i]= 0;
    }

    //Create get and set functions
	public double this[int i]
    {
		get => data[i];
		set => data[i]=value;
	}

    //Simple display function, a list which, for instance, can be written to a table
    public string asList(string sep = "\t")
    {
        string Out = "";
        for (int i = 0; i < size; ++i)
        {
            Out+= this[i];
            if (i<size-1)
                Out+=sep;
        }

        return Out;
    }

    //Pretty display function
    public override string ToString()
    {
        string Out = "";
        for (int i = 0; i < size; ++i)
        {
            Out+= "  |";

            //Round to 5 digits at most, and add some spaces to hopefully keep everything the same length
            Out+=((Abs(this[i])< 100) ? " " : "")+((Abs(this[i])< 10) ? " " : "")+((this[i]>=0) ? " " : "")+string.Format(" {0:N5}",this[i]);

            if (i<size-1)
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
        for (int i = 0; i < size; ++i)
        {

            if (i != size/2)
                for (int j = 0; j<prefix.Length; ++j)
                    Out+= " ";
            else
                Out+=prefix;

            Out+= "|";

            Out+=((Abs(this[i])< 100) ? " " : "")+((Abs(this[i])< 10) ? " " : "")+((this[i]>=0) ? " " : "")+string.Format(" {0:N5}",this[i]);

            if (i<size-1)
                Out+=" |\n";
            else
                Out+= " |";
        }
        return Out;
    }


    //Same as above, but returns a list of strings in order to display matrix multiplication pretty
    public string[] getStrings(string prefix)
    {
        string[] Out = new string[size];
        for (int i = 0; i < size; ++i)
        {
            Out[i]="";
            if (i != size/2)
                for (int j = 0; j<prefix.Length; ++j)
                    Out[i]+= " ";
            else
                Out[i]+=prefix;

            Out[i]+= "|";

            Out[i]+=((Abs(this[i])< 100) ? " " : "")+((Abs(this[i])< 10) ? " " : "")+((this[i]>=-0.001) ? " " : "")+string.Format(" {0:N5}",this[i]);
            Out[i]+=" |";
        }
        return Out;
    }


    public void randomize()
    {
        var generator = new Random();

        for (int i = 0; i < size; ++i)
            this[i]=generator.NextDouble();
    }

    public double[] get_data ()
    {
        return data;
    }


    public double norm ()
    {
        double Norm2 =0;

        for (int i = 0; i < size; ++i)
            Norm2+= this[i]*this[i];
        return Sqrt(Norm2);
    }

    //Essentially the copy constructor which would be used by default in C++
    public vector copy()
    {
        vector Out = new vector(size);
        for(int i=0;i<size;i++)
            Out[i]=this[i];
        return Out;
    }


    //Not as detailed operations as in the matlib class, but enough to get by:
    public static vector operator+(vector a, vector b)
    {
        vector Out = new vector (a.size);
        if (a.size !=b.size )
            throw new ArgumentException("vector size do not match for addition");
        for(int i=0;i<a.size;i++)
            Out[i]=a[i]+b[i];
        return Out;
	}

    public static vector operator-(vector a)
    {
        vector Out = new vector (a.size);

        for(int i=0;i<a.size;i++)
            Out[i]=-a[i];
        return Out;
	}

    public static vector operator- (vector a, vector b)
    {
        vector Out = new vector (a.size);
        if (a.size !=b.size )
            throw new ArgumentException("vector size do not match for addition");
        for(int i=0;i<a.size;i++)
            Out[i]=a[i]-b[i];
        return Out;
	}

    //vector scaling
    public static vector operator*(vector a, double x)
    {
        vector Out = new vector (a.size);
        for(int i=0;i<a.size;i++)
            Out[i]=x*a[i];
        return Out;
    }

    //Either add zeros at the end or delete the last element to get a vector of this length
    public static vector resize(vector a, int n)
    {

        vector Out = new vector (n);
        for(int i=0;i<n;i++)
            Out[i]= i<a.size ? a[i] : 0;
        return Out;
    }

    public static vector operator/(vector a, double x)
    {
        vector Out = new vector (a.size);
        for(int i=0;i<a.size;i++)
            Out[i]=a[i]/x;
        return Out;
    }
    //same idea as in the matrix library given, allow multiplication from both sides
    public static vector operator*(double x, vector a){ return a*x; }

    //Inner product
    public double dot (vector b)
    {
        //Need mathcing matrices
        if (size!=b.size )
            throw new ArgumentException($"vector need to have matching size for inner product to be defined, got {size}  times {b.size} ");

        double Out = 0;
        for (int i=0;i<size;i++)
        {
            Out+=this[i]*b[i];
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

    public bool approx(vector that,double tau=1e-5,double eps = 1e-5)
    {
        if (size != that.size )
            return false;
        for (int i = 0; i < size; ++i)
            if (!approx(that[i],this[i],tau,eps))
                return false;
        return true;
    }


    public bool approx(vector A,vector B,double tau=1e-5,double eps = 1e-5)
    {
        if (A.size != B.size )
            return false;
        for (int i = 0; i < A.size; ++i)
            if (!approx(A[i],B[i],tau,eps))
                return false;
        return true;
    }


}
