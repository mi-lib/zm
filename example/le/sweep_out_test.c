#include <zm/zm_le.h>

void invmat(void)
{
  register int i;
  zMat m1, m2, m;
  zIndex index;
  double marray[] = {
    1, 1, 0, 1,-1,
    0, 1,-1, 3, 0,
   -2,-1, 1, 2, 3,
    1, 0, 2,-1,-1,
   -1,-2,-3,-7, 1
  };

  m1 = zMatCloneArray( marray, 5, 5 );
  m2 = zMatAllocSqr( 5 );
  m  = zMatClone( m1 );
  index = zIndexCreate( 5 );
  zMatIdent( m2 );
  printf( "before sweeped out\n" );
  zMatWrite( m1 );
  zMatWrite( m2 );

  for( i=0; i<zMatRowSize(m1); i++ )
    zSweepOutMat( m1, m2, zIndexElem(index,i), i );
  printf( "after sweeped out\n" );
  zMatWrite( m1 );
  zMatWrite( m2 );

  printf( "m1 * m2 (have to be identity matrix):\n" );
  zMulMatMatDRC( m, m2 );
  zMatWrite( m2 );

  zMatFree( m );
  zMatFree( m1 );
  zMatFree( m2 );
  zIndexFree( index );
}

void leq(void)
{
  register int i;
  zMat m;
  zVec v, x;
  zIndex index;
  double marray[] = {
    1, 1, 2, 1,
    1, 1, 3, 2,
    2,-2, 2,-1,
   -1, 1, 0, 1
  };
  double varray[] = {
   -1, 2,-1, 1
  };

  m = zMatCloneArray( marray, 4, 4 );
  v = zVecCloneArray( varray, 4 );
  x = zVecAlloc( 4 );
  index = zIndexCreate( 4 );
  zIndexWrite( index );
  zMatWrite( m );
  zVecWrite( v );
  zLESolveGauss( m, v, x );

  for( i=0; i<zMatRowSize(m); i++ ){
    zPivoting( m, index, i, zIndexElem(index,i) );
    zSweepOutVec( m, v, zIndexElem(index,i), i );
  }
  printf( "m v =:" );
  zVecWrite( v );
  zIndexWrite( index );
  zVecWrite( x );

  zMatFree( m );
  zVecFree( v );
  zIndexFree( index );
}

int main(void)
{
  invmat();
  zEndl();
  leq();
  return 0;
}
