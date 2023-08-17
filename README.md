ZM - a handy mathematics library.
=================================================================
Copyright (C) Tomomichi Sugihara (Zhidao) since 1998

-----------------------------------------------------------------
## [What is this?]

ZM is a collection of numerical computation and optimization tools
including:

- random number generator (Mersenne twister)
- basic statistics tools
- special functions (Gauss's error, gamma, beta, Bessel)
- complex number (arithmetics, quadratic/cubic equation solvers)
- vectors and matrices of real numbers and complex numbers
- polynomial expression and equation solver
- linear equation solvers (generalized, tridiagonal, Lyapunov)
- matrix decompositions (LU, Choresky, LQ/QR)
- eigenvalue analysis
- finite Fourier series and Fast Fourier Transform
- graph search and rapidly-explored random tree
- kd-tree
- vector sequence and interpolations (linear, Lagrange, Chebyshev,
  spline, Akima, PCHIP, polynomial fitting, NURBS)
- numerical integrator
- optimization solvers (line search, linear programming, linear
  complementary programming, quadratic programming, nonlinear
  programming based on descent/non-descent methods)
- nonlinear equation solvers
- first/second order ordinary differential equation quadratures
  (Runge-Kutta family, embedded Runge-Kutta family, Adams, Gear,
  backward Euler, Butcher-Kuntzmann family, deferred correction,
  Leapflog)
- multiple classification analysis (K-means, GMM)
- growing neural gas algorithm
- nonlinear oscillator
- parser of mathematical expression

ZEDA is required to be installed.

-----------------------------------------------------------------
## [Installation / Uninstallation]

### install

Install ZEDA in advance.

Move to a directly under which you want to install ZM, and run:

   ```
   % git clone https://github.com/zhidao/zm.git
   % cd zm
   ```

Edit **PREFIX** in *config* file if necessary in order to specify
a directory where the header files, the library and some utilities
are installed. (default: ~/usr)

   - header files: $PREFIX/include/zm
   - library file: $PREFIX/lib
   - utilities: $PREFIX/bin

Then, make and install.

   ```
   % make && make install
   ```

### uninstall

Do:

   ```
   % make uninstall
   ```

which removes $PREFIX/lib/libzm.so and $PREFIX/include/zm.

-----------------------------------------------------------------
## [How to use]

When you want to compile your code *test.c*, for example, the following line will work.

   ```
   % gcc `zm-config --cflags` test.c `zm-config -l`
   ```

-----------------------------------------------------------------
## [Contact]

zhidao@ieee.org
