#include <zm/zm_eig.h>

#define TEST 0
#define N    10
int main(void)
{
#if TEST == 1
  int n = 3;
  double a[] = {
    0.0,-6.0, -1.0,
    6.0, 2.0,-16.0,
   -5.0,20.0,-10.0,
  };
  /* eigv = -3.0710, -2.4645+17.6008i, -2.4645-17.6008i
     eigm =
       -0.8326  0.2003-0.1394i   0.2003+0.1394i
       -0.3553 -0.2110-0.6447i  -0.2110+0.6447i
       -0.4248 -0.6930          -0.6930
   */
#elif TEST == 2
  int n = 3;
  double a[] = {
    6.0, 12.0, 19.0,
   -9.0,-20.0,-33.0,
    4.0,  9.0, 15.0,
  };
  /* eigv = -1, 1, 1
     eigm =
       -0.4741   -0.4082   -0.4082
        0.8127    0.8165    0.8165
       -0.3386   -0.4082   -0.4082
   */
#elif TEST == 3
  int n = 3;
  double a[] = {
    3.0, 0.0, 0.0,
    0.0, 2.0,-5.0,
    0.0, 1.0,-2.0,
  };
  /* eigv = 3, i, -i
   */
#elif TEST == 4
  int n = 3;
  double a[] = {
    0, 1, 1,
   -4, 4, 2,
    4,-3,-1,
  };
  /* eigv = 0, 1, 2
   */
#elif TEST == 5
  int n = 3;
  double a[] = {
    1, 0,-1,
   -2, 1, 3,
    2, 1, 2,
  };
  /* eigv = 3, w, w^2
   */
#elif TEST == 6
  int n = 3;
  double a[] = {
    1, 4,-4,
   -1,-3, 2,
    0, 2,-1,
  };
  /* eigv =-3,-1, 1
     eigm =
       -0.8165   0.8165   0.8944
        0.4082  -0.4082   0
       -0.4082  -0.4082   0.4472
   */
#elif TEST == 7
  int n = 3;
  double a[] = {
    2, 1, 0,
    0, 1,-1,
    0, 2, 4,
  };
  /* eigv = 2, 2, 3,
     eigm =
        1  1 -0.4082
        0  0 -0.4082
        0  0  0.8165
   */
#else
  int n = 4;
  double a[] = {
    2, 1, 5, 1,
    1, 3, 7, 0,
    0, 0, 2, 1,
    2, 4, 1, 4,
  };
  /* eigv = 7.03608, 1.38087, 1.29152-2.62195i, 1.29152+2.62195i
   */
#endif
  int i;
  zComplex z[N];
  zMat ma;
  zCMat cma;
  zCVec ve[N], vt;

  ma = zMatCloneArray(a,n,n);
  for( i=0; i<n; i++ )
    ve[i] = zCVecAlloc( n );

  zEigSystem( ma, z, ve, 0 );

  /* ensurance */
  cma = zCMatAlloc( n, n );
  zMat2CMat( ma, cma );
  vt = zCVecAlloc( n );
  for( i=0; i<n; i++ ){
    printf( "eig.#%d val= ", i );
    zComplexWrite( &z[i] );
    printf( "\n" );
    printf( "eig.#%d vec.\n", i );
    zCVecWrite( ve[i] );
    zCMulMatVec( cma, ve[i], vt );
    zComplexRev( &z[i], &z[i] );
    zCVecCatDRC( vt, &z[i], ve[i] );
    printf( " (err) = %g\n", zCVecNorm(vt) );
  }

  zMatFree( ma );
  zCMatFree( cma );
  zCVecFree( vt );
  for( i=0; i<n; i++ )
    zCVecFree( ve[i] );
  return 0;
}
