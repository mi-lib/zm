#include <zm/zm_opt.h>

#define DIM 2

#define TEST 1

double testfunc(zVec x, void *dummy)
{
#if TEST == 0
  return zSqr( zVecElemNC(x,0)-0.4 ) + zSqr( zVecElemNC(x,1)-0.5 );
#else
  /* Goldman & Price function */
  double x1, x2;
  x1 = zVecElemNC(x,0);
  x2 = zVecElemNC(x,1);
  return ( 1 + zSqr(x1+x2+1) * (19-14*x1+3*x1*x1-14*x2+6*x1*x2+3*x2*x2) )
       * ( 30 + zSqr(2*x1-3*x2) * (18-32*x1+12*x1*x1+48*x2-36*x1*x2+27*x2*x2) );
#endif
}

int main(int argc, char *argv[])
{
  zVec ans, min, max;
  double val;

  ans = zVecAlloc( 2 );
  min = zVecCreateList( 2,-2.0,-2.0 );
  max = zVecCreateList( 2, 2.0, 2.0 );
  printf( "iter. num. = %d\n", zOptSolveDIRECT( testfunc, NULL, min, max, 0, zTOL, ans, &val ) );
  zVecPrint( ans );
  printf( "val = %.10g\n", val );
  zVecFreeAO( 3, ans, min, max );
  return 0;
}
