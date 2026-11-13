#include <zm/zm.h>

/* reference values computed by floating-point arithmetics (valid without -ffast-math) */
int ref_isnan(double x){ return x != x; }
int ref_isfinite(double x){ return !ref_isnan( x ) && x - x == 0; }
int ref_isinf(double x){ return ref_isnan( x ) || ref_isfinite( x ) ? 0 : ( x > 0 ? 1 : -1 ); }
int ref_isnormal(double x){ return ref_isfinite( x ) && fabs( x ) >= DBL_MIN; }

/* check all the functions and the macros for a value. */
bool check_ieee(double x)
{
  bool result = true;

  if( zIsInf( x ) != ref_isinf( x ) ){
    eprintf( "zIsInf( %g ) = %d (expected %d)\n", x, zIsInf( x ), ref_isinf( x ) );
    result = false;
  }
  if( zIsNan( x ) != ref_isnan( x ) ){
    eprintf( "zIsNan( %g ) = %d (expected %d)\n", x, zIsNan( x ), ref_isnan( x ) );
    result = false;
  }
  if( zIsFinite( x ) != ref_isfinite( x ) ){
    eprintf( "zIsFinite( %g ) = %d (expected %d)\n", x, zIsFinite( x ), ref_isfinite( x ) );
    result = false;
  }
  if( _zIsInf( x ) != ref_isinf( x ) ){
    eprintf( "_zIsInf( %g ) = %d (expected %d)\n", x, _zIsInf( x ), ref_isinf( x ) );
    result = false;
  }
  if( _zIsNan( x ) != ref_isnan( x ) ){
    eprintf( "_zIsNan( %g ) = %d (expected %d)\n", x, _zIsNan( x ), ref_isnan( x ) );
    result = false;
  }
  if( zIsNormalDouble( x ) != ref_isnormal( x ) ){
    eprintf( "zIsNormalDouble( %g ) = %d (expected %d)\n", x, zIsNormalDouble( x ), ref_isnormal( x ) );
    result = false;
  }
  if( _zIsNormalDouble( x ) != ref_isnormal( x ) ){
    eprintf( "_zIsNormalDouble( %g ) = %d (expected %d)\n", x, _zIsNormalDouble( x ), ref_isnormal( x ) );
    result = false;
  }
  if( _zIsFinite( x ) != ref_isfinite( x ) ){
    eprintf( "_zIsFinite( %g ) = %d (expected %d)\n", x, _zIsFinite( x ), ref_isfinite( x ) );
    result = false;
  }
  return result;
}

/* boundary values of IEEE 754 double-precision floating-point numbers */
bool assert_ieee_boundary(void)
{
  double val[] = {
    0.0, -0.0, 1.0, -1.0,
    DBL_MIN, -DBL_MIN, DBL_MIN/4, -DBL_MIN/4, DBL_MIN*DBL_EPSILON, -DBL_MIN*DBL_EPSILON, /* subnormal numbers */
    DBL_MAX, -DBL_MAX,
    HUGE_VAL, -HUGE_VAL,
  };
  int i;
  bool result = true;

  for( i=0; i<sizeof(val)/sizeof(double); i++ )
    if( !check_ieee( val[i] ) ) result = false;
  return result;
}

/* set the k-th bit from the least significant one of an IEEE 754 double-precision floating-point
 * number, where the 0th to 51st bits are the mantissa part, the 52nd to 62nd bits are the exponent
 * part, and the 63rd bit is the sign bit. */
void set_ieee_bit(__ieee_fp_t *v, int k, bool on)
{
  ubyte *c, mask;

#if __BYTE_ORDER == __BIG_ENDIAN
  c = &v->c[7 - k / 8];
#else
  c = &v->c[k / 8];
#endif
  mask = 1 << ( k % 8 );
  if( on ) *c |= mask; else *c &= ~mask;
}

/* systematic bit patterns: not-a-number with only one bit of the mantissa part set, and numbers
 * with only one bit of the exponent part cleared, for the both signs. */
bool assert_ieee_bits(void)
{
  __ieee_fp_t v;
  int sign, k, j;
  bool result = true;

  for( sign=0; sign<=1; sign++ ){
    for( k=0; k<52; k++ ){ /* not-a-number with a single bit of the mantissa part */
      v.d = 0;
      for( j=52; j<63; j++ ) set_ieee_bit( &v, j, true );
      set_ieee_bit( &v, 63, sign );
      set_ieee_bit( &v, k, true );
      if( !check_ieee( v.d ) ) result = false;
    }
    for( k=52; k<63; k++ ){ /* a finite number with a single bit of the exponent part cleared */
      v.d = 0;
      for( j=52; j<63; j++ ) set_ieee_bit( &v, j, true );
      set_ieee_bit( &v, 63, sign );
      set_ieee_bit( &v, k, false );
      if( !check_ieee( v.d ) ) result = false;
      for( j=0; j<52; j++ ) set_ieee_bit( &v, j, true ); /* with all bits of the mantissa part set */
      if( !check_ieee( v.d ) ) result = false;
    }
    for( k=52; k<63; k++ ){ /* a normal number with a single bit of the exponent part set */
      v.d = 0;
      set_ieee_bit( &v, 63, sign );
      set_ieee_bit( &v, k, true );
      if( !check_ieee( v.d ) ) result = false;
    }
    for( k=0; k<52; k++ ){ /* a subnormal number with a single bit of the mantissa part set */
      v.d = 0;
      set_ieee_bit( &v, 63, sign );
      set_ieee_bit( &v, k, true );
      if( !check_ieee( v.d ) ) result = false;
    }
  }
  return result;
}

/* random bit patterns, half of which have the exponent bits of all one. */
#define N 1000000
bool assert_ieee_random(void)
{
  __ieee_fp_t v;
  int i, j;
  bool result = true;

  for( i=0; i<N; i++ ){
    for( j=0; j<8; j++ ) v.c[j] = zRandI( 0, 0xff );
    if( i % 2 == 0 ){ /* infinity or not-a-number */
      for( j=52; j<63; j++ ) set_ieee_bit( &v, j, true );
      if( i % 4 == 0 ) /* infinity */
        for( j=0; j<52; j++ ) set_ieee_bit( &v, j, false );
    }
    if( !check_ieee( v.d ) ) result = false;
  }
  return result;
}

int main(void)
{
#ifndef __WINDOWS__
  double val;

  zRandInit();
  val = 1.0 / 0.0;
  zAssert( zIsInf, zIsInf( HUGE_VAL ) == 1 &&
                   zIsInf(-HUGE_VAL ) == -1 &&
                   zIsInf( log( 0.0 ) ) == -1 &&
                   zIsInf( val ) == 1 &&
                   zIsInf( NAN ) == 0 &&
                   zIsInf( DBL_MAX ) == 0 &&
                   zIsInf( 0 ) == 0 );
  zAssert( zIsNan, zIsNan( NAN ) &&
                   zIsNan( sqrt( -1 ) ) &&
                   zIsNan( log( -1 ) ) &&
                   !zIsNan( HUGE_VAL ) &&
                   !zIsNan( 0 ) );
  zAssert( zIsFinite, zIsFinite( 1.0 ) &&
                      zIsFinite(-1.0 ) &&
                      zIsFinite( DBL_MAX ) &&
                      !zIsFinite( HUGE_VAL ) &&
                      !zIsFinite(-HUGE_VAL ) &&
                      !zIsFinite( sqrt( -1 ) ) );
  zAssert( zIsInf + zIsNan + zIsFinite + zIsNormalDouble (boundary values), assert_ieee_boundary() );
  zAssert( zIsInf + zIsNan + zIsFinite + zIsNormalDouble (systematic bit patterns), assert_ieee_bits() );
  zAssert( zIsInf + zIsNan + zIsFinite + zIsNormalDouble (random bit patterns), assert_ieee_random() );
#endif /* __WINDOWS__ */
  return EXIT_SUCCESS;
}
