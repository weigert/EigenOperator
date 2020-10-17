# EigenOperator
Generate nD Finite Difference / Lagrange Interpolation / Quadrature Matrix Operators for Eigen, with arbitrary orders and nodes

## What is this repository?

When hand-coding numerical systems, often one needs to construct matrix operators to perform a number of functions:
- Approximation of a derivative at a position based on values at nearby nodes
- Approximation of an integral over a set of nodes
- Approximation of a value at a position based on values at nearby nodes

These operations can be performed using matrix operators.

Using Vandermonde matrices, it is possible to generally derive these nD matrix operators. This library implements a simple algorithm to do just this, and return matrix operators for use with Eigen (C++).

Of course, this is implemented here for regular grids.

## How does it work?
The entire thing is wrapped in a simple namespace that you can use to generate the matrix operators you need for your system.

## ToDo
General boundary conditions and non-periodic boundary conditions.

## License
MIT License
