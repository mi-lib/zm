#include <zm/zm.h>

#define DIM 7
#define TOL 1.0e-7

void assert_peqsolve(void)
{
  double a[DIM];
  zPex pex;
  zCVec ans;
  zComplex c;
  register int i;
  bool ret = true;

  for( i=0; i<DIM; i++ ) a[i] = zRandF(-5,5);
  pex = zPexCloneArray( a, DIM-1 );
  ans = zCVecAlloc( DIM );

  zPexDKA( pex, ans, zTOL, 0 );
  for( i=0; i<DIM-1; i++ ){
    zPexCVal( pex, zCVecElemNC(ans,i), &c );
    if( !zComplexIsTol( &c, zTOL ) ) ret = false;
  }
  zAssert( zPexDKA, ret );

  zPexBH( pex, ans, zTOL, 0 );
  for( i=0; i<DIM-1; i++ ){
    zPexCVal( pex, zCVecElemNC(ans,i), &c );
    if( !zComplexIsTol( &c, TOL ) ) ret = false;
  }
  zAssert( zPexBH, ret );
  zCVecFree( ans );
}

#define NUM 7

void assert_exp(void)
{
  zPex p;
  zVec fact;
  zCVec ans;
  register int i;
  bool ret = true;

  fact = zVecAlloc( NUM );
  ans = zCVecAlloc( NUM );
  for( i=0; i<NUM; i++ )
    zVecElemNC(fact,i) = zRandF(-10,10);
  p = zPexExp( fact );
  zPexDKA( p, ans, zTOL, 0 );
  for( i=0; i<NUM; i++ ){
    if( !zComplexIsReal( zCVecElemNC(ans,i), zTOL ) ) ret = false;
    if( !zVecValIsIncluded( fact, zCVecElemNC(ans,i)->re, zTOL ) ) ret = false;
  }
  zPexFree( p );
  zVecFree( fact );
  zCVecFree( ans );
  zAssert( zPexExp, ret );
}

#define NUM_REAL 3
#define NUM_IMAG 2

void assert_cexp(void)
{
  zPex p;
  zCVec fact, ans;
  zComplex c;
  register int i;
  bool ret = true;

  zRandInit();
  fact = zCVecAlloc( NUM_REAL + NUM_IMAG * 2 );
  for( i=0; i<NUM_REAL; i++ )
    zComplexCreate( zCVecElemNC(fact,i), zRandF(-10,10), 0 );
  for( i=0; i<NUM_IMAG; i++ ){
    zComplexCreate( zCVecElemNC(fact,2*i+NUM_REAL), zRandF(-10,10), zRandF(-10,10) );
    zComplexConj( zCVecElemNC(fact,2*i+NUM_REAL), zCVecElemNC(fact,2*i+NUM_REAL+1) );
  }
  p = zPexCExp( fact );

  ans = zCVecAlloc( zCVecSizeNC(fact) );
  zPexDKA( p, ans, TOL, 0 );
  zCVecTouchup( ans );

  for( i=0; i<zCVecSizeNC(fact); i++ ){
    zPexCVal( p, zCVecElemNC(fact,i), &c );
    if( !zComplexIsTol( &c, TOL ) ) ret = false;
  }
  zPexFree( p );
  zCVecFree( fact );
  zCVecFree( ans );
  zAssert( zPexCExp, ret );
}

int main(void)
{
  zRandInit();
  assert_peqsolve();
  assert_exp();
  assert_cexp();
  return EXIT_SUCCESS;
}
