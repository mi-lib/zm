#include <zm/zm_eig.h>

#define TEST 4
int main(void)
{
#if TEST == 1
  double a[] = {
    2, 3,-2,
   -3,14,-7,
   -5,19,-9,
  };
  int n = 3, m = 3;
  /* ans =
     | 0.11 -0.97  0.21 || 27.02 0    0    || -0.21  0.88 -0.43 |
     | 0.59  0.11 -0.80 ||  0    2.82 0    ||  0.96 -0.10  0.26 |
     | 0.80  0.22  0.56 ||  0    0    0.16 ||  0.19  0.46  0.87 |
   */
#elif TEST == 2
  double a[] = {
    9, 4,
    6, 8,
    2, 7,
  };
  int n = 3, m = 2;
  /* ans =
     | -0.6105  0.7174  0.3355 || 14.9359 0      || -0.6925 -0.7214 |
     | -0.6646 -0.2336 -0.7098 ||  0      5.1883 ||  0.7214 -0.6925 |
     | -0.4308 -0.6563  0.6194 ||  0      0      |
   */
#elif TEST == 3
  double a[] = {
    1, 2,
    3, 4,
    5, 6,
    7, 8,
  };
  int n = 4, m = 2;
  /* ans =
 | 0.1525  0.8226 -0.3945   -0.3800 ||
 | 0.3499  0.4214  0.2428    0.8007 || 14.2691 0      || 0.6414 -0.7672 |
 | 0.5474  0.0201  0.6979   -0.4614 ||  0      0.6268 || 0.7672  0.6414 |
 | 0.7448 -0.3812 -0.5462    0.0407 ||  0      0      |
   */
#else
  double a[] = {
    1, 2, 3, 4, 5,
    5, 6, 7, 8, 9,
    9,10,11,12,13,
  };
  int n = 3, m = 5;
#endif
  zMat ma, u, v, l, tmp1, tmp2;
  zVec sv;
  int i, rank;

  ma = zMatCloneArray( a, n, m );
  u = zMatAllocSqr( n );
  v = zMatAlloc( n, m );
  sv = zVecAlloc( n );

  zMatWrite( ma );
  rank = zSVD( ma, sv, u, v );
  zVecWrite( sv );
  zMatWrite( u );
  zMatWrite( v );

  printf( ">>ensurance\n" );
  l = zMatAlloc( n, rank );
  tmp1 = zMatAlloc( n, rank );
  tmp2 = zMatAlloc( n, m );
  for( i=0; i<rank; i++ )
    zMatSetElem( l, i, i, zVecElem(sv,i) );
  zMulMatMat( u, l, tmp1 );
  zMulMatMat( tmp1, v, tmp2 );
  zMatSubDRC( tmp2, ma );
  printf( "|| ULV - A || = %g\n", zMatNorm(tmp2) );

  zMatFree( ma );
  zMatFree( u );
  zMatFree( v );
  zMatFree( l );
  zMatFree( tmp1 );
  zMatFree( tmp2 );
  zVecFree( sv );
  return 0;
}
