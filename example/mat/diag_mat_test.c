#include <zm/zm_mat.h>

int main(void)
{
  zMat m;
  zVec diag;

  diag = zVecAlloc( 5 );
  zVecLinSpace( diag, 0, 100 );

  m = zMatAllocSqr( zVecSize(diag) );
  printf( "%dx%d identity matrix.\n", zVecSize(diag), zVecSize(diag) );
  zMatIdent( m );
  zMatWrite( m );
  printf( "%dx%d diagonal matrix.\n", zVecSize(diag), zVecSize(diag) );
  zVecWrite( diag );
  zMatDiag( m, diag );
  zMatWrite( m );

  zMatFree( m );
  zVecFree( diag );
  return 0;
}
