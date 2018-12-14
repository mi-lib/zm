#include <zm/zm_eig.h>

#define ROW 20
#define COL 8

int main(void)
{
  zMat ma, u, v, l, tmp1, tmp2;
  zVec sv;
  int i, rank;

  ma = zMatAlloc( ROW, COL );
  u = zMatAllocSqr( ROW );
  v = zMatAlloc( ROW, COL );
  sv = zVecAlloc( ROW );

  zRandInit();
  zMatRandUniform( ma, -10, 10 );
  zMatWrite( ma );
  rank = zSVD( ma, sv, u, v );
  zVecWrite( sv );
  zMatWrite( u );
  zMatWrite( v );

  printf( ">>ensurance\n" );
  l = zMatAlloc( ROW, rank );
  tmp1 = zMatAlloc( ROW, rank );
  tmp2 = zMatAlloc( ROW, COL );
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
