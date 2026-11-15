#include <zm/zm.h>
#include <float.h>

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
  int i, n;
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
    zPexFree( q );
    zPexFree( r );
    zPexFree( g );
  }
  zAssert( zPexDiv + zPexMul + zPexAdd, result );
  zPexFree( p );
  zPexFree( f );
}

/* check if the residuals of a polynomial at the roots are not greater than tol times the scales of
 * the rounding errors of the values, sum_k |a_k| |z|^k. */
#define TOL_PEX_EQ ( 100 * DBL_EPSILON ) /* tentative */
bool check_pex_roots(const zPex a, const zCVec ans, double tol)
{
  int i, k;
  double az, scale;
  zComplex c;
  bool result = true;

  for( i=0; i<zPexDim(a); i++ ){
    zPexCVal( a, zCVecElemNC(ans,i), &c );
    az = zComplexAbs( zCVecElemNC(ans,i) );
    for( scale=0, k=zPexDim(a); k>=0; k-- ) scale = scale * az + fabs( zPexCoeff(a,k) );
    if( !( zComplexAbs( &c ) <= tol * scale ) ){
      eprintf( "root #%d: relative residual = %g\n", i, zComplexAbs( &c ) / scale );
      result = false;
    }
  }
  return result;
}

/* check if the roots are the same with the expected ones relatively. */
bool check_pex_roots_expected(const zCVec ans, const double *re, int n, double tol)
{
  int i, j;
  double d, dmin;
  zComplex c;
  bool result = true;

  for( i=0; i<n; i++ ){
    for( dmin=HUGE_VAL, j=0; j<n; j++ ){
      zComplexCreate( &c, zCVecElemNC(ans,j)->re - re[i], zCVecElemNC(ans,j)->im );
      if( ( d = zComplexAbs( &c ) ) < dmin ) dmin = d;
    }
    if( !( dmin <= tol * fabs( re[i] ) ) ){
      eprintf( "root %g not found (error = %g)\n", re[i], dmin );
      result = false;
    }
  }
  return result;
}

void assert_peqsolve(void)
{
  double a[DIM];
  zPex pex;
  zCVec ans;
  int i;
  bool ret_dka = true, ret_bh = true;

  for( i=0; i<DIM; i++ ) a[i] = zRandF(-5,5);
  pex = zPexCloneArray( a, DIM-1 );
  ans = zCVecAlloc( DIM );
  if( !zPexDKA( pex, ans, 0, 0 ) || !check_pex_roots( pex, ans, TOL_PEX_EQ ) ) ret_dka = false;
  if( !zPexBH( pex, ans, 0, 0 ) || !check_pex_roots( pex, ans, TOL_PEX_EQ ) ) ret_bh = false;
  zPexFree( pex );
  zCVecFree( ans );
  zAssert( zPexDKA, ret_dka );
  zAssert( zPexBH, ret_bh );
}

/* polynomials with roots of various scales, and a polynomial with a small leading coefficient. */
void assert_peqsolve_scale(void)
{
  const double scale[] = { 1.0e-60, 1.0e-13, 1.0e-6, 1.0, 1.0e3, 1.0e60 }; /* scale^4 has to be representable */
  double re[4];
  zPex pex;
  zCVec ans;
  int i, k;
  bool ret_dka = true, ret_bh = true;

  pex = zPexAlloc( 4 );
  ans = zCVecAlloc( 4 );
  for( k=0; k<(int)( sizeof(scale)/sizeof(double) ); k++ ){
    for( i=0; i<4; i++ ) re[i] = ( i + 1 ) * scale[k];
    /* (x-s)(x-2s)(x-3s)(x-4s) */
    zPexSetCoeffList( pex, 24*pow(scale[k],4), -50*pow(scale[k],3), 35*scale[k]*scale[k], -10*scale[k], 1.0 );
    if( !zPexDKA( pex, ans, 0, 0 ) || !check_pex_roots_expected( ans, re, 4, 1.0e-12 ) ) ret_dka = false;
    if( !zPexBH( pex, ans, 0, 0 ) || !check_pex_roots_expected( ans, re, 4, 1.0e-12 ) ) ret_bh = false;
  }
  /* (x+2)^4, for which the initial points of Durand-Kerner-Aberth's method coincide */
  for( i=0; i<4; i++ ) re[i] = -2;
  zPexSetCoeffList( pex, 16.0, 32.0, 24.0, 8.0, 1.0 );
  if( !zPexDKA( pex, ans, 0, 0 ) || !check_pex_roots_expected( ans, re, 4, 1.0e-12 ) ) ret_dka = false;
  if( !zPexBH( pex, ans, 0, 0 ) || !check_pex_roots_expected( ans, re, 4, 1.0e-3 ) ) ret_bh = false; /* a quadruple root */
  /* a small leading coefficient, which makes a large root */
  zPexSetCoeffList( pex, -1.0, 2.0, -3.0, 4.0, 1.0e-3 );
  if( !zPexDKA( pex, ans, 0, 0 ) || !check_pex_roots( pex, ans, TOL_PEX_EQ ) ) ret_dka = false;
  if( !zPexBH( pex, ans, 0, 0 ) || !check_pex_roots( pex, ans, TOL_PEX_EQ ) ) ret_bh = false;
  zPexFree( pex );
  zCVecFree( ans );
  zAssert( zPexDKA (various scales), ret_dka );
  zAssert( zPexBH (various scales), ret_bh );
}

#define NUM 7

void assert_exp(void)
{
  zPex p;
  zVec fact;
  zCVec ans;
  int i;
  bool ret = true;

  fact = zVecAlloc( NUM );
  ans = zCVecAlloc( NUM );
  for( i=0; i<NUM; i++ )
    zVecElemNC(fact,i) = zRandF(-10,10);
  p = zPexExp( fact );
  zPexDKA( p, ans, 0, 0 );
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
  int i;
  bool ret = true;

  fact = zCVecAlloc( NUM_REAL + NUM_IMAG * 2 );
  for( i=0; i<NUM_REAL; i++ )
    zComplexCreate( zCVecElemNC(fact,i), zRandF(-10,10), 0 );
  for( i=0; i<NUM_IMAG; i++ ){
    zComplexCreate( zCVecElemNC(fact,2*i+NUM_REAL), zRandF(-10,10), zRandF(-10,10) );
    zComplexConj( zCVecElemNC(fact,2*i+NUM_REAL), zCVecElemNC(fact,2*i+NUM_REAL+1) );
  }
  p = zPexCExp( fact );

  ans = zCVecAlloc( zCVecSizeNC(fact) );
  zPexDKA( p, ans, 0, 0 );
  for( i=0; i<zCVecSizeNC(fact); i++ ){
    zPexCVal( p, zCVecElemNC(fact,i), &c );
    if( !zComplexIsTol( &c, TOL ) ) ret = false;
    if( !zComplexValIsIncluded( zCVecBufNC(ans), zCVecSizeNC(ans), zCVecElemNC(fact,i), TOL ) ) ret = false;
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
  assert_peqsolve_scale();
  assert_exp();
  assert_cexp();
  return EXIT_SUCCESS;
}
