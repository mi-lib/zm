ZM - a handy mathematics library.
Copyright (C) 1998 Tomomichi Sugihara (Zhidao)

-----------------------------------------------------------------
[What is this?]

ZM is a collection of numerical computation and optimization
tools including:

 - random number generator (Mersenne twister)
 - statistics
 - special functions (Gauss's error, gamma, beta, Bessel)
 - vectors and matrices of real numbers and complex numbers
 - complex number (arithmetics, quadratic/cubic equation solvers)
 - polynomial expression and equation solver
 - linear equation solvers (generalized, tridiagonal, Lyapunov)
 - matrix decompositions (LU, Choresky, LQ/QR)
 - eigenvalue analysis
 - Fast Fourier Transform
 - graph search and rapidly-explored random tree.
 - kd-tree
 - vector sequence and interpolations (linear, Lagrange, spline
   Akima, polynomial fitting, NURBS)
 - numerical integrator
 - optimization tools (line search, linear programming, linear
   complementary programming, quadratic programming, nonlinear
   programming based on descent/non-descent methods)
 - nonlinear equation solvers
 - first/second order ordinary differential equation quadratures
   (Runge-Kutta family, embedded Runge-Kutta family, Adams, Gear,
   backward Euler, Butcher-Kuntzmann family, deferred correction,
   Leapflog)
 - multiple classification analysis (clustering)
 - nonlinear oscillator
 - parser of mathematical expression

ZEDA is required.

-----------------------------------------------------------------
[Installation / Uninstallation]

<install>
0. Install ZEDA in advance.

1. Unpack the distributed archive where you want.

% zcat zm-X.Y.Z.tgz | tar xvf
or, if you use GNU tar,
% tar xzvf zm-X.Y.Z.tgz

X.Y.Z is for the revision.

2. Enter the directory.

% cd zm-X.Y.Z

3. Edit config file if necessary.
  PREFIX   directory where the library is installed.
           ~/usr as a default. In this case, header files
           and library are installed under ~/usr/include
           and ~/usr/lib, respectively.

4. Make it.

% make

5. Install it.

% make install

Or,

% cp -a lib/libzm.so $PREFIX/lib/
% cp -a include/zm $PREFIX/include/
% cp -a bin/* $PREFIX/bin/

<uninstall>
Delete $PREFIX/lib/libzm.so and $PREFIX/include/zm.

-----------------------------------------------------------------
[How to use]

You may need to set your PATH and LD_LIBRARY_PATH environment
variables. This is done by:
 export PATH=$PATH:$PREFIX/bin
 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PREFIX/lib
if your working shell is Bourne shell (bash, zsh, etc.), or by:
 set path = ( $path $PREFIX/bin )
 setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH:$PREFIX/lib
if your working shell is C shell (csh, tcsh, etc.).

When you want to compile your code test.c, for example, the following
line will work.

% gcc `zm-config -L` `zm-config -I` test.c `zm-config -l`

-----------------------------------------------------------------
[Contact]

zhidao@ieee.org
