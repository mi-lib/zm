#include <zm/zm.h>

#define DIM 7
#define TOL ( 1.0e-6 )

double pex_val1(zPex p, double arg)
{
  int i, n;
  double val = 0, b;

  if( arg == 0 ) return zPexCoeff( p, 0 );
  n = zPexDim(p);
  if( fabs(arg) < 1.0 )
    for( i=n, b=pow(arg,n); i>=0; i--, b/=arg )
      val += zPexCoeff(p,i) * b;
  else
    for( i=0, b=1.0; i<=n; i++, b*=arg )
      val += zPexCoeff(p,i) * b;
  return val;
}

double pex_val2(zPex p, double arg)
{
  int i, n;
  double val;

  if( arg == 0 ) return zPexCoeff( p, 0 );
  n = zPexDim(p);
  for( val=zPexCoeff(p,n), i=n-1; i>=0; i-- )
    val = val * arg + zPexCoeff(p,i);
  return val;
}

double pex_val3(zPex p, double arg)
{
  int i, n;
  double val, b, *term;
  zIndex index;

  if( arg == 0 ) return zPexCoeff( p, 0 );
  n = zPexDim(p);
  if( !( term = zAlloc( double, n+1 ) ) ) return NAN;
  if( !( index = zIndexAlloc( n+1 ) ) ){
    zIndexFree( index );
    return NAN;
  }
  zIndexOrder( index, 0 );
  for( i=0, b=1.0; i<=n; i++, b*=arg )
    term[i] = zPexCoeff(p,i) * b;
  zDataSortAbsIndex( term, n+1, index );
  for( val=0, i=n; i>=0; i-- )
    val += term[zIndexElemNC(index,i)];
  zIndexFree( index );
  free( term );
  return val;
}

void assert_pex_val(void)
{
  zPex p;
  double x;
  register int i, n;
  bool result = true;

  p = zPexCreateList( 3, zRandF(-1,1), zRandF(-1,1), zRandF(-1,1), zRandF(-1,1) );
  n = 100;
  for( i=0; i<=n; i++ ){
    x = 5.0*i/n - 2.0;
    if( !zEqual( zPexVal(p,x), pex_val3(p,x), zTOL ) ) result = false;
  }
  zPexFree( p );
  zAssert( zPexVal, result );
}

#define N 100

void assert_pex_cval(void)
{
  int i;
  double t;
  zPex p;
  zComplex arg, c;
  bool result = true;

  p = zPexAlloc( N );
  zPexSetCoeff( p, N, 1.0 );
  zPexSetCoeff( p, 0,-1.0 );
  for( i=0; i<=N; i++ ){
    t = 2 * zPI * i / N;
    zComplexCreatePolar( &arg, 1.0, t );
    zPexCVal( p, &arg, &c );
    if( !zComplexIsTiny( &c ) ) result = false;
  }
  zPexFree( p );
  zAssert( zPexCVal, result );
}

void assert_pex_div(void)
{
  zPex p, f, q, r, g;
  int i;
  bool result = true;

  p = zPexAlloc( 5 );
  f = zPexAlloc( 3 );
  for( i=0; i<N; i++ ){
    zVecRandUniform( p, -10, 10 );
    zVecRandUniform( f, -10, 10 );
    zPexDiv( p, f, &q, &r );
    g = zPexMul( f, q );
    zPexAddDRC( g, r );
    if( !zBoolStr( zPexEqual( p, g, zTOL ) ) ) result = false;
  }
  zAssert( zPexDiv + zPexMul + zPexAdd, result );
  zPexFree( p );
  zPexFree( f );
  zPexFree( q );
  zPexFree( r );
}

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
    if( !zComplexIsTol( &c, TOL ) ){
      zComplexPrint( &c ); zEndl(); fflush( stdout );
      ret = false;
    }
  }
  zAssert( zPexDKA, ret );

  zPexBH( pex, ans, zTOL, 0 );
  for( i=0; i<DIM-1; i++ ){
    zPexCVal( pex, zCVecElemNC(ans,i), &c );
    if( !zComplexIsTol( &c, TOL ) ){
      zComplexPrint( &c ); zEndl(); fflush( stdout );
      ret = false;
    }
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
    if( !zComplexIsReal( zCVecElemNC(ans,i), TOL ) ) ret = false;
    if( !zVecValIsIncluded( fact, zCVecElemNC(ans,i)->re, TOL ) ) ret = false;
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
  zPexDKA( p, ans, zTOL, 0 );
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
  assert_pex_val();
  assert_pex_cval();
  assert_pex_div();
  assert_peqsolve();
  assert_exp();
  assert_cexp();
  return EXIT_SUCCESS;
}
