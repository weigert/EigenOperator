# EigenOperator
Single-Header n-Dimensional Differentiation / Integration / Interpolation Matrix Operators for Eigen. Arbitrary Approximation Orders, Arbitrary Support Nodes.

## What is this repository?

When hand-coding numerical systems, often one needs to construct matrix operators to perform a number of functions:
- Approximation of a derivative at a position based on values at nearby nodes
- Approximation of an integral over a set of nodes
- Approximation of a value at a position based on values at nearby nodes

This allows for full-discretization of continuous equations, and their solution as a linear algebra problem. This is generally useful for e.g. the solution of ordinary or partial differential equations.

This repository will algorithmically generate these matrix operators as a simple Eigen::SparseMatrix. It also works for arbitrary dimensions of data, by flattening the nD-Array to a 1D system and constructing the matrix appropriately.

The order of approximation and the choice of support points can be arbitrary. Construction of the matrix operators is very simple and intuitive.

### How it Works
By inverting the vandermonde matrix for a given set of support points, we can compute the weights required to interpolate a value from the support points (i.e. solve the linear system).

Similarly, by inverting the taylor coefficient matrix for a given set of k points, we receive a matrix of weights required to approximate all derivatives around the origin up to the k'th derivative order, to maximum accuracy.

The quadrature approximation weights have not yet been implemented.

The set of weights are then converted into matrix operators to yield the interpolation and differentiation approximations at all points based on the chosen discretization scheme.

This is all implemented for regular grids with periodic boundary conditions, for simplicity reasons. The non-periodic boundary conditions can be handled by properly composing multiple periodic matrices.

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
