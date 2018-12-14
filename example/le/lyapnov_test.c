#include <zm/zm_le.h>

#define TEST 0

int main(void)
{
#if TEST == 1
  double ar[] = { 2, -1, 5, 1 };
  double br[] = { 47, 8, -13, -11 };
  int n = 2;
  double xr[] = { 3, 7, 0, -2 }; /* true answer */
#elif TEST == 2
  double ar[] = { 1, -3, 2, 1 };
  double br[] = { -6, -5, -5, 14 };
  int n = 2;
  double xr[] = { 1, -2, -2, 1 }; /* true answer */
#elif TEST == 3
  double ar[] = { 1, 0, 2, 2, 1,-1, -1, 3, 1 };
  double br[] = { 11, 4, 8, 3, 10, 9, 2, 1, -2 };
  int n = 3;
  double xr[] = { 3, 1, -1, 1, 2, 3, 0, -1, 1 }; /* true answer */
#else
  double ar[] = { 1, -2, 3, 1 };
  double br[] = { 7, 1, 3, -8 };
  int n = 2;
  double xr[] = { -1, 1, 2, -1 }; /* true answer */
#endif
  zMat a, b, x;

  a = zMatCloneArray( ar, n, n );
  b = zMatCloneArray( br, n, n );
  x = zMatAllocSqr( n );
  zLyapnovSolve( a, b, x );
  printf( "A = " ); zMatWrite( a );
  printf( "B = " ); zMatWrite( b );
  zMatTouchup( x );
  zMatWrite( x );
  zMatCopyArray( xr, n, n, x );
  printf( "X(true) = " ); zMatWrite( x );
  return 0;
}
