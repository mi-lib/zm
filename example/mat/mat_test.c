#include <zm/zm_mat.h>

#define R 3
#define C 4

int main(void)
{
  zMat m, mc;
  zVec diag, row, col;
  double val[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

  m = zMatAlloc( R, C );
  mc = zMatAlloc( R, C );
  diag = zVecCreateList( 3, -100.0, -200.0, -300.0 );

  printf( "clear\n" );
  zMatClear( m );
  zMatWrite( m );
  printf( "identity\n" );
  zMatIdentNC( m );
  zMatWrite( m );
  printf( "diag\n" );
  zMatDiagNC( m, diag );
  zMatWrite( m );

  printf( "clone\n" );
  mc = zMatClone( m );
  zMatWrite( mc );
  printf( "copy array\n" );
  zMatCopyArray( val, R, C, mc );
  zMatWrite( mc );

  printf( "abstract of a row and col\n" );
  row = zVecAlloc( C );
  col = zVecAlloc( R );
  zMatGetRow( mc, 1, row );
  zVecWrite( row );
  zMatGetCol( mc, 1, col );
  zVecWrite( col );
  printf( "replace of a row and col\n" );
  zMatSetRow( m, 0, row );
  zMatWrite( m );
  zMatSetCol( m, 2, col );
  zMatWrite( m );
  printf( "trace\n" );
  printf( "%.10f\n", zMatTrNC( m ) );
  printf( "norm\n" );
  printf( "%.10f\n", zMatNorm( m ) );

  zMatFree( m );
  return 0;
}
