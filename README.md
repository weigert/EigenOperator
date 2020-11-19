# EigenOperator
Single-Header n-Dimensional Differentiation / Integration / Interpolation Matrix Operators for Eigen. Arbitrary Approximation Orders, Arbitrary Support Nodes, Arbitrary Boundary Conditions.

Use this to discretize, linearize and solve your continuous equations systems!

## What is this repository?

When hand-coding numerical systems, often one needs to construct matrix operators to perform a number of functions:
- Approximation of a derivative at a position based on values at nearby nodes
- Approximation of an integral over a set of nodes
- Approximation of a value at a position based on values at nearby nodes

This allows for full-discretization of continuous equations, and their solution as a linear algebra problem. This is generally useful for e.g. the solution of ordinary or partial differential equations and a large variety of other systems!

This repository will algorithmically generate these matrix operators as a simple Eigen::SparseMatrix. It also works for arbitrary dimensions of data, by flattening the nD-Array to a 1D system and constructing the matrix appropriately.

The order of approximation and the choice of support points can be arbitrary. Construction of the matrix operators is very simple and intuitive.

## Usage

Simply include the headerfile "operator.h" in your Eigen project.

      #include "operator.h"
      
      //...
      
EigenOperator is wrapped in the namespace `op` that contains all the functions needed to generate the matrix operators you need for your system.
      
### Matrix Operators

Matrix operators are constructed by generating a set of weights based on the support points and the operation required. The weights are then associated with grid elements and the matrix operator is constructed.

#### Finite Differences

Finite difference operators require an offset to the neighboring points from which the differential will be approximated and a differentiation order. Note that the differentiation order must be strictly smaller than the number of support points. The approximation order is mathematically defined based on the differentiation order and number of supports.

      //Finite Differences (Examples)
      std::vector<double> weights = op::FD({support_node_offsets}, differentiation_order); 
      
      //2nd order central differences, 2nd order approximation
      std::vector<double> weights = op::FD({-1, 0, 1}, 2);
      
      //2nd order derivative, central differences, 3rd order approximation
      std::vector<double> weights = op::FD({-2, -1, 0, 1, 2}, 2);
      
      //1st order derivative, forward differences, 2nd order approximation 
      std::vector<double> weights = op::FD({0, 1, 2}, 1);
      
#### Lagrange Interpolation

Lagrange interpolation operators generate an approximation of a continuous value from a known set of neighboring values.
      
      //Lagrange Interpolation (Examples)
      std::vector<double> weights = op::LI({support_node_offsets});
      
      //Average between two neighboring nodes
      std::vector<double> weights = op::LI({-0.5, 0.5});
      
      //Interpolate value from 3 forward values
      std::vector<double> weights = op::LI({1.0, 2.0, 3.0});
      
#### Quadrature / Integration
      
To-Do

#### Operator Construction

Once a set of weights has been generated for a given operator, the operator is constructed by association with other vector elements:
      
      //Previously: typedef SparseMatrix<double, RowMajor> sparse;

      vec dim = {X, Y, Z, W, ...}; //Size of the grid (any dimension number)
      initializer_list<vec> points = {...}; //Support point positions (relative to centroid where we evaluate)
      sparse MATRIX_OPERATOR = op::make(points, weights, dim);
    
The vector "dim" gives the size of the grid in n-Dimensions. points represents a vector of relative positions (offsets) to the support points, for which the weights were computed. The individual vectors have the dimension of the grid, so you can e.g. construct a diffusion operator in only a single spatial direction.

We finally pass the weights and the matrix is filled appropriately.

### Boundary Conditions

By default, the operators are constructed with periodic boundary conditions. To generate alternative boundary conditions, matrix operators can be combined using masking matrices.

A masking matrix is constructed using a lambda function of an n-dimensional position vector p, which returns a boolean, setting the value to 0 or 1.

	sparse bcmask = op::mask(dim, [&](vec p){
		return some_test_of_p;
	});
      
This can be used to e.g. use lower approximation orders at the boundary.

### How it Works

By inverting the vandermonde matrix for a given set of support points, we can compute the weights required to interpolate a value from the support points (i.e. solve the linear system).

Similarly, by inverting the taylor coefficient matrix for a given set of k points, we receive a matrix of weights required to approximate all derivatives around the origin up to the k'th derivative order, to maximum accuracy.

The quadrature approximation weights have not yet been implemented.

The set of weights are then converted into matrix operators to yield the interpolation and differentiation approximations at all points based on the chosen discretization scheme.

This is all implemented for regular grids with periodic boundary conditions, for simplicity reasons. The non-periodic boundary conditions can be handled by properly composing multiple periodic matrices.

## ToDo

The finite difference and lagrange interpolator operators currently only support 1D weight evaluation, i.e. you cannot lagrange interpolate using a 2D polynomial currently. This will probably be expanded in the future.

## License
MIT License
