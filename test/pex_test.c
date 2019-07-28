#include <zm/zm.h>

#define DIM 7
#define TOL 1.0e-10

void assert_peqsolve(void)
{
  double a[DIM];
  zPex pex;
  zComplex ans[DIM], c;
  register int i;
  bool ret = true;

  for( i=0; i<DIM; i++ ) a[i] = zRandF(-5,5);
  pex = zPexCloneArray( a, DIM-1 );

  zPexDKA( pex, ans, zTOL, 0 );
  for( i=0; i<DIM-1; i++ ){
    zPexCVal( pex, &ans[i], &c );
    if( !zComplexIsTol( &c, TOL ) ) ret = false;
  }
  zAssert( zPexDKA, ret );

  zPexBH( pex, ans, zTOL, 0 );
  for( i=0; i<DIM-1; i++ ){
    zPexCVal( pex, &ans[i], &c );
    if( !zComplexIsTol( &c, TOL ) ) ret = false;
  }
  zAssert( zPexBH, ret );
}

int main(void)
{
  zRandInit();
  assert_peqsolve();
  return EXIT_SUCCESS;
}
