#include <zm/zm_mat_eig.h>

#define N 100

bool assert_mat_eig_one(double a[], int n)
{
  int i;
  zCVec eigval, eigvec, vt;
  zMat ma;
  zCMat cma, eigbase;
  zComplex z;
  bool result = true;
  const double tol = 1.0e-7;

  ma = zMatCloneArray( a, n, n );
  cma = zCMatAllocSqr( n );
  eigval = zCVecAlloc( n );
  eigbase = zCMatAllocSqr( n );
  eigvec = zCVecAlloc( n );
  vt = zCVecAlloc( n );
  zMatEig( ma, eigval, eigbase, 0 );
  /* confirmation */
  zMatToCMat( ma, cma );
  for( i=0; i<n; i++ ){
    zCMatGetCol( eigbase, i, eigvec );
    zCMulMatVec( cma, eigvec, vt );
    zComplexRev( zCVecElemNC(eigval,i), &z );
    zCVecCatDRC( vt, &z, eigvec );
    if( !zCVecIsTol( vt, tol ) ){
      eprintf( " error norm = %g\n", zCVecNorm(vt) );
      result = false;
    }
  }
  zMatFree( ma );
  zCMatFree( cma );
  zCMatFree( eigbase );
  zCVecFree( eigval );
  zCVecFree( eigvec );
  zCVecFree( vt );
  return result;
}

void assert_mat_eig_random(void)
{
  zMat ma;
  const int n = 20;
  int i;
  bool result = true;

  ma = zMatAllocSqr( n );
  for( i=0; i<N; i++ ){
    zMatRandUniform( ma, -10, 10 );
    if( !assert_mat_eig_one( zMatBuf(ma), n ) ) result = false;
  }
  zMatFree( ma );
  zAssert( zMatEig (random), result );
}

int main(void)
{
  zRandInit();
  assert_mat_eig_random();
  return 0;
}
