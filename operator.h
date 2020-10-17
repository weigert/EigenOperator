/*

EigenOperator
Author: Nicholas McDonald

//

There are three main operator types:

  - Interpolators   (i.e. lagrange polynomials)
  - Differentiators (i.e. finite differences)
  - Integrators     (i.e. numerical quadrature) - not implemented at the moment

As they are linearizations, they can be put into general matrix form after
the individual operator weights have been computed.

Below are weight calculating functions for all three types of operators, as
well as one general function for placing them into a sparse matrix system.

This assumes a uniform grid structure!

*/

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

typedef Triplet<double> triplet;
typedef SparseMatrix<double, RowMajor> sparse;
typedef Array<int, 1, Dynamic> vec;

/*
===============================
  Overloaded Helper Operators
===============================
*/

vec operator%(const vec& rhs, const vec& lhs){
  vec n(lhs.size());
  for(int i = 0; i < rhs.size(); i++)
    n[i] = rhs[i]%lhs[i];
  return n;
}

ostream& operator<<(ostream& out, const vec data){
  out<<"[";
  for(unsigned int i = 0; i < data.size()-1; i++)
    out<<data[i]<<", ";
  out<<data[data.size()-1]<<"]";
  return out;
}

template<typename T>
ostream& operator<<(ostream& out, const vector<T> data){
  for(unsigned int i = 0; i < data.size(); i++)
    out<<data[i]<<" ";
  return out;
}

namespace op{

/*
=============================
    nD Space Linearization
=============================
*/

vec itop(int i, vec& dim){
  vec n(dim.size());
  for(unsigned int j = 0; j < dim.size(); j++){
    int F = 1;
    for(unsigned int k = dim.size()-1; k > j; k--)
      F *= dim[j];
    n[j] = (i/F)%dim[j];
  }
  return n;
}

int ptoi(vec& pos, vec& dim){
  vec w(dim.size());
  for(unsigned int i = 0; i < dim.size(); i++){
    int F = 1;
    for(unsigned int j = dim.size()-1; j > i; j--)
      F *= dim[j];
    w[i] = F;
  }
  return (pos*w).sum();
}

/*
=============================
    Operator Construction
=============================
*/

sparse make(initializer_list<vec> _p, vector<double>& w, vec& dim){

  //Assert that dim > 0
  dim.cwiseMin(0);

  unsigned int SIZE = dim.prod(); //Compute total Vector Size
  sparse M(SIZE, SIZE);

  vector<triplet> list; //Triplet List for Matrix Construction
  vector<vec> p = _p;   //Position Vector Extracted

  for(unsigned int i = 0; i < SIZE; i++){                     //Iterate over Positions
    for(unsigned int j = 0; j < p.size(); j++){     //Iterate over Nodes

      //Compute Position of Shifted Element, add to list with weights
      vec shift(dim.size());
      shift = (itop(i, dim) + dim + p[j]);
      shift = shift % dim;
      std::cout<<shift<<" "<<ptoi(shift, dim)<<std::endl;
      
      list.push_back(triplet(i, ptoi(shift, dim), w[j]));

    }
  }

  //Construct Opeator Matrix and Return
  M.setFromTriplets(list.begin(), list.end());
  std::cout<<M<<std::endl;

  return M;
}

/*
=============================
      Weight Computation
=============================
*/

//Factorial Evaluator
double fac(int k){
  int j = 1;  //Factorial
  for(int i = 1; i <= k; i++)
    j*= i;
  return j;
}

//Taylor Coefficients
double taylor(double x, double a, int n){
  return (pow(x-a, n) / fac(n));
}

/*
  Lagrange Interpolator:
    Setup the Vandermonde Matrix for the linear equation of the polynomial.
    Inversion and evaluation at point = 0 yields the linear combination of
    evaluated values f, that gives us point a0 (the intercept), which is
    our interpolated value for arbitrary polynomial degree.
*/

vector<double> LI(initializer_list<double> shift){
  vector<double> weights = shift;

  //Vandermonde Matrix
  int N = shift.size();
  MatrixXd V(N, N);
  for(int p = 0; p < N; p++)
    for(int d = 0; d < N; d++)
      V(p, d) = pow(weights[p], d);

  //Invert the Matrix
  V = V.inverse();

  //Extract the 0th Row
  for(int i = 0; i < N; i++)
    weights[i] = V(0, i);

  return weights;
}

/*
  Finite Difference Approximator:
    Compute the Taylor Matrix for all points expanded around zero
    Invert the matrix and select the appropriate row for the weights!
*/

vector<double> FD(initializer_list<double> points, unsigned int order){
  vector<double> weights = points;

  //Check for Consistency
  if(order >= points.size()){
    cout<<"Order must be strictly smaller than the number of support points."<<endl;
    return weights;
  }

  //Taylor Matrix
  int N = points.size();
  MatrixXd T(N, N);
  for(int p = 0; p < N; p++)
    for(int d = 0; d < N; d++)
      T(p, d) = taylor(weights[p], 0, d);

  //Invert the Matrix
  T = T.inverse();

  //Extract the Order'th Row
  for(int i = 0; i < N; i++)
    weights[i] = T(order, i);

  return weights;
}

/*
  To Do: Numerical Quadrature Approximator
*/

/*
  General Sparse Operator Matrix:

  Quadrature, Finite Differences and Lagrange Interpolation are all linearized forms of the mathematical integration, differentiation and interpolation operators. They are represented my matrices operating on a state vector.

  Each operator has a set of weights which are computed for the relevant operator type. They matrices are then constructed from the weight sets.

  We specify the nodes over which we wish to compute our approximation. Then the weights are computed and the matrix operator is constructed.
*/

/*
Required data:

-> Nodes, their position in space (relative)
-> The set of weights associated with each node
-> The overall size of the domain
*/



/*
================================================================================
                      EXAMPLE OPERATORS FOR CONVENIENCE
================================================================================
*/

/*
  These individual operators can be tuned to any accuracy you like.
*/

/*

//Interpolate the surface values linearly on both faces in a direction!
SparseMatrix<double,RowMajor> FV_FLUX(vec2 dim){
  vector<double> w = LI({-0.5, 0.5});
  return OP({0, 1}, w, dim) -  OP({0, -1}, w, dim);
  //vector<double> w = FD({-1, 0, 1}, 1);
  //return OP({-1, 0, 1}, w, dim);
}

//1st Order Accurate Forward Finite Differences, Differential at Boundary for
//both surfaces in a direction summed up!
SparseMatrix<double,RowMajor> FV_DIFFUSION(vec2 dim){
  vector<double> w1 = FD({0, 1}, 1);
  vector<double> w2 = FD({0, -1}, 1);
  return OP({0, 1}, w1, dim) -  OP({0, -1}, w2, dim);
}

//4th Order Accurate Centered Differences, Approximation of Nth order Derivative
SparseMatrix<double,RowMajor> CFD(vec2 dim, int n){
  vector<double> w = FD({-2, -1, 0, 1, 2}, n);
  return OP({-2, -1, 0, 1, 2}, w, dim);
}

//Forward Finite Differences
SparseMatrix<double,RowMajor> FFD(vec2 dim, int n){
  vector<double> w = FD({0, 1}, n);
  return OP({0, 1}, w, dim);
}

//Forward Finite Differences
SparseMatrix<double,RowMajor> BFD(vec2 dim, int n){
  vector<double> w = FD({0, 1}, n);
  return OP({-1, 0}, w, dim);
}

//Forward Finite Differences
SparseMatrix<double,RowMajor> SFD(vec2 dim, int n){
  vector<double> w = FD({-1, 0, 1}, n);
  return OP({-1, 0, 1}, w, dim);
}

*/

}; //End of Namespace
