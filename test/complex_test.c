#include <zm/zm.h>
#include <float.h>

#define NUM 100

/* the backward error of an answer x of a polynomial equation coeff[n] x^n + ... + coeff[0] = 0, namely,
 * the residual relative to its scale sum_k |coeff[k]||x|^k (the residual itself if the scale is zero).
 * for |x| > 1, both are divided by x^n and evaluated as polynomials of 1/x so as not to overflow. */
double pe_backward_error(const double coeff[], int n, zComplex *x)
{
  zComplex r, rw, w;
  double aw, s;
  int k;

  if( zComplexAbs( x ) > 1 ){
    zComplexInv( x, &w );
    aw = zComplexAbs( &w );
    zComplexCreate( &r, coeff[0], 0 );
    for( s=fabs(coeff[0]), k=1; k<=n; k++ ){
      zComplexCMul( &w, &r, &rw );
      zComplexCreate( &r, rw.re + coeff[k], rw.im );
      s = s * aw + fabs( coeff[k] );
    }
  } else{
    aw = zComplexAbs( x );
    zComplexCreate( &r, coeff[n], 0 );
    for( s=fabs(coeff[n]), k=n-1; k>=0; k-- ){
      zComplexCMul( x, &r, &rw );
      zComplexCreate( &r, rw.re + coeff[k], rw.im );
      s = s * aw + fabs( coeff[k] );
    }
  }
  return s > 0 ? zComplexAbs( &r ) / s : zComplexAbs( &r );
}

#define TOL_PE_SOLVE ( 100*DBL_EPSILON )

bool check_pe_answers(const char *name, const double coeff[], int n, zComplex ans[])
{
  int i;
  double e;
  bool result = true;

  for( i=0; i<n; i++ )
    if( !( ( e = pe_backward_error( coeff, n, &ans[i] ) ) <= TOL_PE_SOLVE ) ){
      eprintf( "%s: backward error %g of %g%+gi for coefficients", name, e, ans[i].re, ans[i].im );
      for( i=n; i>=0; i-- ) eprintf( " %.17g", coeff[i] );
      eprintf( "\n" );
      result = false;
      break;
    }
  return result;
}

bool check_qesolve(double a, double b, double c)
{
  double coeff[] = { c, b, a };
  zComplex ans[2];

  zQESolve( a, b, c, ans );
  return check_pe_answers( "zQESolve", coeff, 2, ans );
}

bool check_cesolve(double a, double b, double c, double d)
{
  double coeff[] = { d, c, b, a };
  zComplex ans[3];

  zCESolve( a, b, c, d, ans );
  if( ans[0].im != 0 ){
    eprintf( "zCESolve: the first answer is not real for (%.17g, %.17g, %.17g, %.17g)\n", a, b, c, d );
    return false;
  }
  return check_pe_answers( "zCESolve", coeff, 3, ans );
}

/* random coefficients, which include a small leading coefficient, coefficients making a nearly double
 * answer, and coefficients of various scales. */
bool assert_qesolve(void)
{
  int i;
  double a, b, c;
  bool result = true;

  for( i=0; i<NUM; i++ ){
    a = zRandF(-3,3); if( a == 0 ) a = 1;
    b = zRandF(-3,3);
    c = zRandF(-3,3);
    if( !check_qesolve( a, b, c ) ) result = false;
    if( !check_qesolve( ( a > 0 ? 1 : -1 ) * pow( 10, zRandF(-12,-4) ), b, c ) ) result = false;
    if( !check_qesolve( a, b, b*b/(4*a)*( 1 + zRandF(-1.0e-8,1.0e-8) ) ) ) result = false;
    if( !check_qesolve( a*pow( 10, zRandF(-150,150) ), b*pow( 10, zRandF(-150,150) ), c*pow( 10, zRandF(-150,150) ) ) ) result = false;
  }
  return result;
}

/* equations with zero, multiple and imaginary answers */
bool assert_qesolve_special(void)
{
  struct{
    double a, b, c;
  } *eq, special[] = {
    { 2,  0,  0 }, /* 2 x^2 */
    { 1, -1,  0 }, /* x ( x - 1 ) */
    { 1, -2,  1 }, /* ( x - 1 )^2 */
    { 1,  0,  1 }, /* ( x - i )( x + i ) */
    { 1,  0, -1 }, /* ( x - 1 )( x + 1 ) */
    { 0,  0,  0 },
  };
  bool result = true;

  for( eq=special; eq->a!=0; eq++ )
    if( !check_qesolve( eq->a, eq->b, eq->c ) ) result = false;
  return result;
}

/* random coefficients, which include
 * - a small leading coefficient, for which the answers except the largest one suffered from cancellation,
 * - c ~ b^2/(3a), for which the depressed cubic equation has p ~ 0 and Cardano's formula suffered from
 *   cancellation, and
 * - coefficients of various scales. */
bool assert_cesolve(void)
{
  int i;
  double a, b, c, d;
  bool result = true;

  for( i=0; i<NUM; i++ ){
    a = zRandF(-3,3); if( a == 0 ) a = 1;
    b = zRandF(-3,3);
    c = zRandF(-3,3);
    d = zRandF(-3,3);
    if( !check_cesolve( a, b, c, d ) ) result = false;
    if( !check_cesolve( ( a > 0 ? 1 : -1 ) * pow( 10, zRandF(-12,-4) ), b, c, d ) ) result = false;
    if( !check_cesolve( a, b, b*b/(3*a)*( 1 + zRandF(-1.0e-6,1.0e-6) ), d ) ) result = false;
    if( !check_cesolve( a*pow( 10, zRandF(-8,8) ), b*pow( 10, zRandF(-8,8) ), c*pow( 10, zRandF(-8,8) ), d*pow( 10, zRandF(-8,8) ) ) ) result = false;
    /* the ratio of the magnitudes of the answers is within the limit about 1e+150 (see zCESolve()) */
    if( !check_cesolve( a*pow( 10, zRandF(-30,30) ), b*pow( 10, zRandF(-30,30) ), c*pow( 10, zRandF(-30,30) ), d*pow( 10, zRandF(-30,30) ) ) ) result = false;
  }
  return result;
}

/* equations with zero, multiple and imaginary answers */
bool assert_cesolve_special(void)
{
  struct{
    double a, b, c, d;
  } *eq, special[] = {
    { 1,  0,  0,  0 }, /* x^3 */
    { 1, -1,  0,  0 }, /* x^2 ( x - 1 ) */
    { 1,  0, -1,  0 }, /* x ( x - 1 )( x + 1 ) */
    { 1,  0,  1,  0 }, /* x ( x - i )( x + i ) */
    { 1, -3,  3, -1 }, /* ( x - 1 )^3 */
    { 1, -3,  0,  4 }, /* ( x - 2 )^2 ( x + 1 ) */
    { 1,  0,  0, -1 }, /* x^3 - 1 */
    { 1, -6, 11, -6 }, /* ( x - 1 )( x - 2 )( x - 3 ) */
    { 0,  0,  0,  0 },
  };
  bool result = true;

  for( eq=special; eq->a!=0; eq++ )
    if( !check_cesolve( eq->a, eq->b, eq->c, eq->d ) ) result = false;
  return result;
}

bool assert_isincluded(void)
{
  int i;
  zComplex c[NUM];

  for( i=0; i<NUM; i++ )
    zComplexCreate( &c[i], zRandF(-10,10), zRandF(-10,10) );
  for( i=0; i<NUM; i++ )
    if( !zComplexValIsIncluded( c, NUM, &c[i], zTOL ) ) return false;
  return true;
}

bool assert_div_drc(void)
{
  zComplex c1, c2, c3;
  double k;
  int i;

  for( i=0; i<NUM; i++ ){
    zComplexCreate( &c1, zRandF(-10,10), zRandF(-10,10) );
    k = zRandF( 1, 10 );
    zComplexDiv( &c1, k, &c2 );
    c3 = c1;
    zComplexDivDRC( &c3, k );
    if( !zComplexEqual( &c2, &c3, zTOL ) ) return false;
    c3 = c1;
    _zComplexDivDRC( &c3, k );
    if( !zComplexEqual( &c2, &c3, zTOL ) ) return false;
  }
  return true;
}

/* check if q c2 = c1 with the relative residual less than 10 eps. */
bool check_cdiv(const zComplex *c1, const zComplex *c2, const zComplex *q)
{
  zComplex r;

  zComplexCMul( q, c2, &r );
  zComplexSubDRC( &r, c1 );
  return zComplexAbs( &r ) <= 10 * DBL_EPSILON * zComplexAbs( c1 );
}

/* division of complex numbers with very large or small absolute values */
bool assert_cdiv(void)
{
  zComplex c1, c2, q, inv, one;
  double scale[] = { 1, 1.0e-7, 1.0e-100, 1.0e-200, 1.0e-300, 1.0e100, 1.0e200, 1.0e300 };
  int i, j;

  zComplexCreate( &one, 1, 0 );
  for( i=0; i<NUM; i++ )
    for( j=0; j<sizeof(scale)/sizeof(double); j++ ){
      zComplexCreate( &c1, zRandF(-10,10), zRandF(-10,10) );
      zComplexCreate( &c2, zRandF(-10,10)*scale[j], zRandF(-10,10)*scale[j] );
      if( !zComplexCDiv( &c1, &c2, &q ) || !check_cdiv( &c1, &c2, &q ) ){
        eprintf( "(%g%+gi)/(%g%+gi) = %g%+gi\n", c1.re, c1.im, c2.re, c2.im, q.re, q.im );
        return false;
      }
      if( zComplexCDivDRC( &c1, &c2 ) != &c1 || !zComplexEqual( &c1, &q, 0 ) ) return false;
      if( fabs( c2.re ) < 1.0e-150 || fabs( c2.re ) > 1.0e150 ) continue; /* 1/c2 might overflow */
      if( !zComplexInv( &c2, &inv ) || !check_cdiv( &one, &c2, &inv ) ) return false;
    }
  /* purely real and purely imaginary divisors */
  zComplexCreate( &c1, 3, 4 );
  zComplexCreate( &c2, 0, 2 );
  if( !zComplexCDiv( &c1, &c2, &q ) || q.re != 2 || q.im != -1.5 ) return false;
  zComplexCreate( &c2, 2, 0 );
  if( !zComplexCDiv( &c1, &c2, &q ) || q.re != 1.5 || q.im != 2 ) return false;
  return true;
}

/* division by real numbers with very small or large absolute values, and normalization */
bool assert_div_normalize_tol(void)
{
  zComplex c, q, n;
  double k[] = { 3.0, 1.0e-13, 4.0e-309, 1.0e308 };
  int j;

  for( j=0; j<sizeof(k)/sizeof(double); j++ ){
    zComplexCreate( &c, 1.0e-300 * zRandF(-1,1), 1.0e-300 * zRandF(-1,1) );
    if( !zComplexDiv( &c, k[j], &q ) ) return false;
    if( fabs( q.re - c.re/k[j] ) > DBL_EPSILON * fabs( c.re/k[j] ) ||
        fabs( q.im - c.im/k[j] ) > DBL_EPSILON * fabs( c.im/k[j] ) ) return false;
  }
  if( zComplexDiv( &c, 0, &q ) ) return false;
  zComplexCreate( &c, 3.0e-14, 4.0e-14 ); /* absolute value 5e-14 */
  if( zComplexNormalize( &c, 0, &n ) || n.re != 0 || n.im != 0 ) return false;
  if( !zComplexNormalize( &c, 1.0e-20, &n ) || !zIsTol( zComplexAbs( &n ) - 1, 10*DBL_EPSILON ) ) return false;
  return true;
}

/* division by zero */
bool assert_cdiv_zero(void)
{
  zComplex c1, c2, q;

  zComplexCreate( &c1, 1, 2 );
  zComplexCreate( &c2, 0, 0 );
  q = c1;
  return !zComplexCDiv( &c1, &c2, &q ) && !zComplexInv( &c2, &q ) &&
         !zComplexCDivDRC( &q, &c2 ) && zComplexEqual( &q, &c1, 0 );
}

/* multiplications where the result is put into an operand */
bool assert_cmul_inplace(void)
{
  zComplex c1, c2, c3, c4;
  int i;

  for( i=0; i<NUM; i++ ){
    zComplexCreate( &c1, zRandF(-10,10), zRandF(-10,10) );
    zComplexCreate( &c2, zRandF(-10,10), zRandF(-10,10) );
    zComplexCMul( &c1, &c2, &c3 );
    c4 = c1;
    zComplexCMul( &c4, &c2, &c4 );
    if( !zComplexEqual( &c3, &c4, zTOL ) ) return false;
    c4 = c2;
    zComplexCMul( &c1, &c4, &c4 );
    if( !zComplexEqual( &c3, &c4, zTOL ) ) return false;
    zComplexCMulConj( &c1, &c2, &c3 );
    c4 = c1;
    zComplexCMulConj( &c4, &c2, &c4 );
    if( !zComplexEqual( &c3, &c4, zTOL ) ) return false;
  }
  return true;
}

/* absolute values of complex numbers, the squared absolute values of which underflow or overflow. */
bool assert_abs_range(void)
{
  const int scale[] = { -1022, -600, 0, 600, 1022 };
  zComplex c0, c;
  double v, ref;
  int i, k;
  bool result = true;

  for( i=0; i<1000; i++ ){
    /* the real part, which dominates the absolute value, is not less than 0.5 */
    _zComplexCreate( &c0, zRandF(0.5,1) * ( zRandI(0,1) ? 1 : -1 ), zRandF(-1,1) );
    ref = sqrt( c0.re*c0.re + c0.im*c0.im );
    for( k=0; k<(int)( sizeof(scale)/sizeof(int) ); k++ ){
      _zComplexCreate( &c, ldexp( c0.re, scale[k] ), ldexp( c0.im, scale[k] ) );
      v = ldexp( zComplexAbs( &c ), -scale[k] );
      if( !( fabs( v - ref ) <= 4 * DBL_EPSILON * ref ) ) result = false;
    }
  }
  _zComplexCreate( &c, ldexp( 1.0, -1074 ), 0 ); /* the smallest subnormal number */
  if( zComplexAbs( &c ) != ldexp( 1.0, -1074 ) ) result = false;
  _zComplexCreate( &c, DBL_MAX, 0 );
  if( zComplexAbs( &c ) != DBL_MAX ) result = false;
  _zComplexCreate( &c, 0, 0 );
  if( zComplexAbs( &c ) != 0 ) result = false;
  _zComplexCreate( &c, HUGE_VAL, 1 );
  if( zComplexAbs( &c ) != HUGE_VAL ) result = false;
  return result;
}

/* relative judgements of real numbers and co-conjugate numbers, which do not depend on the scale. */
bool assert_isreal_coconj(void)
{
  const double scale[] = { 1.0e-300, 1.0e-100, 1.0, 1.0e100, 1.0e300 };
  zComplex c1, c2;
  int k;

  for( k=0; k<(int)( sizeof(scale)/sizeof(double) ); k++ ){
    _zComplexCreate( &c1, scale[k], 1.0e-13 * scale[k] );
    if( !zComplexIsReal( &c1, zTOL ) ) return false;
    _zComplexCreate( &c1, scale[k], 1.0e-11 * scale[k] );
    if( zComplexIsReal( &c1, zTOL ) ) return false;
    _zComplexCreate( &c1, scale[k], 2 * scale[k] );
    _zComplexCreate( &c2, scale[k], -2 * ( 1 + 1.0e-13 ) * scale[k] );
    if( !zComplexCoconj( &c1, &c2, zTOL ) ) return false;
    _zComplexCreate( &c2, scale[k], -2 * ( 1 + 1.0e-11 ) * scale[k] );
    if( zComplexCoconj( &c1, &c2, zTOL ) ) return false;
  }
  zComplexZero( &c1 );
  zComplexZero( &c2 );
  return zComplexIsReal( &c1, zTOL ) && zComplexCoconj( &c1, &c2, zTOL );
}

int main(int argc, char *argv[])
{
  zRandInit();
  zAssert( zComplexDivDRC + _zComplexDivDRC, assert_div_drc() );
  zAssert( zComplexCMul + zComplexCMulConj (in-place), assert_cmul_inplace() );
  zAssert( zComplexCDiv + zComplexCDivDRC + zComplexInv, assert_cdiv() );
  zAssert( zComplexDiv + zComplexNormalize (tolerance), assert_div_normalize_tol() );
  zAssert( zComplexCDiv + zComplexCDivDRC + zComplexInv (division by zero), assert_cdiv_zero() );
  zAssert( zComplexValIsIncluded, assert_isincluded() );
  zAssert( zComplexAbs (underflow and overflow), assert_abs_range() );
  zAssert( zComplexIsReal + zComplexCoconj (relative judgements), assert_isreal_coconj() );
  zAssert( zQESolve (random coefficients of various scales), assert_qesolve() );
  zAssert( zQESolve (zero, multiple and imaginary answers), assert_qesolve_special() );
  zAssert( zCESolve (random coefficients of various scales), assert_cesolve() );
  zAssert( zCESolve (zero, multiple and imaginary answers), assert_cesolve_special() );
  return 0;
}
