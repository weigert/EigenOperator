# EigenOperator
Single-Header n-Dimensional Differentiation / Integration / Interpolation Matrix Operators for Eigen. Arbitrary Approximation Orders, Arbitrary Support Nodes.

## What is this repository?

When hand-coding numerical systems, often one needs to construct matrix operators to perform a number of functions:
- Approximation of a derivative at a position based on values at nearby nodes
- Approximation of an integral over a set of nodes
- Approximation of a value at a position based on values at nearby nodes

These operations can be performed using matrix operators.

Using Vandermonde matrices, it is possible to generally derive these nD matrix operators. This library implements a simple algorithm to do just this, and return matrix operators for use with Eigen (C++).

Of course, this is implemented here for regular grids.

## Usage

The entire thing is wrapped in a simple namespace "op" that you can use to generate the matrix operators you need for your system. For a specific operator, we generate a set of weights using the required support nodes. We also specify the derivative / approximation order:

      //Finite Differences (Examples)
      
      //2nd order central differences, 2nd order approximation
      std::vector<double> weights = op::FD({-1, 0, 1}, 2);
      
      //2nd order derivative, central differences, 3rd order approximation
      std::vector<double> weights = op::FD({-2, -1, 0, 1, 2}, 2);
      
      //1st order derivative, forward differences, 2nd order approximation 
      std::vector<double> weights = op::FD({0, 1, 2}, 1);
      
      
      //Lagrange Interpolation (Examples)
      
      //Average between two neighboring nodes
      std::vector<double> weights = op::LI({-0.5, 0.5});
      
      //Interpolate value from 3 forward values
      std::vector<double> weights = op::LI({1.0, 2.0, 3.0});
      
We then take our weight set and construct our operator:

      //Previously: typedef SparseMatrix<double, RowMajor> sparse;

      vec dim = {X, Y, Z, W, ...}; //Size of the grid (any dimension number)
      initializer_list<vec> points = {...}; //Support point positions (relative to centroid where we evaluate)
      sparse MATRIX_OPERATOR = op::make(points, weights, dim);
    
The vector "dim" gives the size of the grid in n-Dimensions. points represents a vector of relative positions (offsets) to the support points, for which the weights were computed. The individual vectors have the dimension of the grid, so you can e.g. construct a diffusion operator in only a single spatial direction.

We finally pass the weights and the matrix is filled appropriately.

## ToDo
General boundary conditions and non-periodic boundary conditions.

The finite difference and lagrange interpolator operators currently only support 1D weight evaluation, i.e. you cannot lagrange interpolate using a 2D polynomial currently. This will probably be expanded in the future.

## License
MIT License
