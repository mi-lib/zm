#include <zm/zm.h>

#define TEST_VEC_SIZE    10
#define TEST_VEC_BUFSIZ 100

void assert_raw_vec_get_put(void)
{
  double test_vec1[TEST_VEC_BUFSIZ];
  double test_vec2[TEST_VEC_BUFSIZ];
  double test_vec3[TEST_VEC_BUFSIZ];

  zRawVecRandUniform( test_vec1, TEST_VEC_SIZE, -10, 10 );
  zRawVecRandUniform( test_vec2, TEST_VEC_SIZE, -10, 10 );
  zRawVecGet( test_vec1, 2, test_vec3, 3 );
  zAssert( zRawVecGet, memcmp( test_vec1+2, test_vec3, sizeof(double)*3 ) == 0 );
  zRawVecPut( test_vec2, 2, test_vec3, 3 );
  zAssert( zRawVecPut, memcmp( test_vec2+2, test_vec3, sizeof(double)*3 ) == 0 );
  zRawVecCopy( test_vec1, test_vec2, TEST_VEC_SIZE );
  zRawVecSwap( test_vec2, 3, 6 );
  zAssert( zRawVecSwap, test_vec1[3] == test_vec2[6] && test_vec1[6] == test_vec2[3] );
}

void assert_raw_vec_misc(void)
{
  double test_vec1[TEST_VEC_BUFSIZ];
  double test_vec2[TEST_VEC_BUFSIZ];
  double val1, val2, dval;
  int i;
  bool result;

  val1 = zRandF( -10, 10 );
  val2 = zRandF( -10, 10 );
  zRawVecSetAll( test_vec1, TEST_VEC_SIZE, val1 );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( test_vec1[i] != val1 ) result = false;
  zAssert( zRawVecSetAll, result );
  zRawVecLinSpace( test_vec1, TEST_VEC_SIZE, val1, val2 );
  dval = ( val2 - val1 ) / ( TEST_VEC_SIZE - 1 );
  for( result=true, i=1; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i] - test_vec1[i-1] - dval ) ) result = false;
  zAssert( zRawVecLinSpace, result );
  zRawVecShift( test_vec1, TEST_VEC_SIZE, dval, test_vec2 );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec2[i] - test_vec1[i] - dval ) ) result = false;
  zAssert( zRawVecShift, result );
  zRawVecCopy( test_vec1, test_vec2, TEST_VEC_SIZE );
  zRawVecShiftDRC( test_vec2, TEST_VEC_SIZE, dval );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec2[i] - test_vec1[i] - dval ) ) result = false;
  zAssert( zRawVecShiftDRC, result );
}

void assert_raw_vec_equal(void)
{
  double *v1, *v2;
  const int size = 10;
  int i;
  bool result1, result2;

  v1 = zAlloc( double, size );
  v2 = zAlloc( double, size );
  if( !v1 || !v2 ){
    ZALLOCERROR();
    return;
  }
  zRawVecRandUniform( v1, size, -10, 10 );
  zRawVecCopy( v1, v2, size );
  for( i=0; i<size; i++ )
    *( v2 + i ) += zTOL * 0.5;
  result1 = zRawVecEqual( v1, v2, size, zTOL );
  for( i=0; i<size; i++ )
    *( v2 + i ) += zTOL * 10;
  result2 = zRawVecEqual( v1, v2, size, zTOL );
  free( v1 );
  free( v2 );
  zAssert( zRawVecEqual (positive case), result1 );
  zAssert( zRawVecEqual (negative case), !result2 );
}

void assert_raw_vec_match(void)
{
  double *v1, *v2;
  const int size = 10;
  bool result;

  v1 = zAlloc( double, size );
  v2 = zAlloc( double, size );
  if( !v1 || !v2 ){
    ZALLOCERROR();
    return;
  }
  zRawVecRandUniform( v1, size, -10, 10 );
  zRawVecCopy( v1, v2, size );
  result = zRawVecMatch( v1, v2, size );
  free( v1 );
  free( v2 );
  zAssert( zRawVecCopy + zRawVecMatch, result );
}

void assert_raw_vec_arith(void)
{
  double test_vec1[TEST_VEC_BUFSIZ];
  double test_vec2[TEST_VEC_BUFSIZ];
  double test_vec3[TEST_VEC_BUFSIZ];
  double k;
  int i;
  bool result;

  zRawVecRandUniform( test_vec1, TEST_VEC_SIZE, -10, 10 );
  zRawVecRandUniform( test_vec2, TEST_VEC_SIZE, -10, 10 );
  k = zRandF( -10, 10 );

  zRawVecAdd( test_vec1, test_vec2, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]+test_vec2[i]-test_vec3[i] ) ) result = false;
  zAssert( zRawVecAdd, result );
  zRawVecSub( test_vec1, test_vec2, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]-test_vec2[i]-test_vec3[i] ) ) result = false;
  zAssert( zRawVecSub, result );
  zRawVecRev( test_vec1, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]+test_vec3[i] ) ) result = false;
  zAssert( zRawVecRev, result );
  zRawVecMul( test_vec1, k, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]*k-test_vec3[i] ) ) result = false;
  zAssert( zRawVecMul, result );
  zRawVecDiv( test_vec1, k, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]/k-test_vec3[i] ) ) result = false;
  zAssert( zRawVecDiv, result );
  zRawVecAmp( test_vec1, test_vec2, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]*test_vec2[i]-test_vec3[i] ) ) result = false;
  zAssert( zRawVecAmp, result );
  zRawVecDem( test_vec1, test_vec2, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]/test_vec2[i]-test_vec3[i] ) ) result = false;
  zAssert( zRawVecDem, result );
  zRawVecCat( test_vec1, k, test_vec2, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]+k*test_vec2[i]-test_vec3[i] ) ) result = false;
  zAssert( zRawVecCat, result );
}

void assert_raw_vec_innerprod(void)
{
  const double v1[]  = { 3,-2, 5 };
  const double v2[]  = {-6, 1,-4 };
  const double v3[]  = { 5, 3,-6 };
  const double vec[] = { 2,-1, 3 };
  const int size = 3;

  zAssert( zRawVecInnerProd,
    zEqual( zRawVecInnerProd( v1, vec, size ), 23, 0 ) &&
    zEqual( zRawVecInnerProd( v2, vec, size ),-25, 0 ) &&
    zEqual( zRawVecInnerProd( v3, vec, size ),-11, 0 ) );
}

void assert_raw_vec_normalize(void)
{
  double test_vec1[TEST_VEC_BUFSIZ];
  double test_vec2[TEST_VEC_BUFSIZ];
  double norm;
  int i;
  bool result;

  zRawVecRandUniform( test_vec1, TEST_VEC_SIZE, -10, 10 );
  norm = zRawVecNorm( test_vec1, TEST_VEC_SIZE );
  zRawVecNormalize( test_vec1, TEST_VEC_SIZE, 0, test_vec2 );
  zAssert( zRawVecNormalize, zIsTiny( zRawVecNorm(test_vec2,TEST_VEC_SIZE) - 1 ) );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]/test_vec2[i] - norm ) ) result = false;
  zAssert( zRawVecNorm, result );
}

/* normalization of tiny vectors with tolerances */
void assert_raw_vec_normalize_tol(void)
{
  double v1[TEST_VEC_BUFSIZ], v2[TEST_VEC_BUFSIZ];
  int i;
  bool result_default = true, result_tol = true, result_nan;

  zRawVecRandUniform( v1, TEST_VEC_SIZE, -1, 1 );
  zRawVecMulDRC( v1, 1.0e-13 / zRawVecNorm( v1, TEST_VEC_SIZE ), TEST_VEC_SIZE ); /* norm = 1e-13 */
  /* rejected with the default tolerance, and the result is the zero vector */
  zRawVecRandUniform( v2, TEST_VEC_SIZE, -1, 1 );
  if( zRawVecNormalize( v1, TEST_VEC_SIZE, 0, v2 ) ) result_default = false;
  for( i=0; i<TEST_VEC_SIZE; i++ )
    if( v2[i] != 0 ) result_default = false;
  zRawVecRandUniform( v2, TEST_VEC_SIZE, -1, 1 );
  if( zRawVecNormalize( v1, TEST_VEC_SIZE, -1, v2 ) ) result_default = false; /* negative tolerance */
  /* normalized with a smaller tolerance */
  if( !zRawVecNormalize( v1, TEST_VEC_SIZE, 1.0e-20, v2 ) ||
      !zIsTol( zRawVecNorm( v2, TEST_VEC_SIZE ) - 1, 10*DBL_EPSILON ) ) result_tol = false;
  /* rejected with a larger tolerance */
  if( zRawVecNormalize( v1, TEST_VEC_SIZE, 1.0e-10, v2 ) ) result_tol = false;
  /* not-a-number */
  v1[0] = NAN;
  result_nan = zRawVecNormalize( v1, TEST_VEC_SIZE, 0, v2 ) == NULL;
  zAssert( zRawVecNormalize (default tolerance), result_default );
  zAssert( zRawVecNormalize (specified tolerance), result_tol );
  zAssert( zRawVecNormalize (not-a-number), result_nan );
}

/* division by very small or large values */
void assert_raw_vec_div_range(void)
{
  double v1[TEST_VEC_BUFSIZ], v2[TEST_VEC_BUFSIZ];
  double k[] = { 3.0, 1.0e-13, 4.0e-309, 1.0e-320, 1.0e308, -1.0e308 }; /* including subnormal numbers */
  int i, j;
  bool result = true, result_drc = true;

  for( j=0; j<sizeof(k)/sizeof(double); j++ ){
    zRawVecRandUniform( v1, TEST_VEC_SIZE, -1, 1 );
    if( k[j] < 1.0e-300 ) zRawVecMulDRC( v1, 1.0e-300, TEST_VEC_SIZE ); /* avoid overflow of quotients */
    zRawVecDiv( v1, k[j], v2, TEST_VEC_SIZE );
    for( i=0; i<TEST_VEC_SIZE; i++ ) /* relative error not greater than one ulp */
      if( !zIsFinite( v2[i] ) || fabs( v2[i] - v1[i]/k[j] ) > DBL_EPSILON * fabs( v1[i]/k[j] ) ) result = false;
    zRawVecDivDRC( v1, k[j], TEST_VEC_SIZE );
    for( i=0; i<TEST_VEC_SIZE; i++ )
      if( v1[i] != v2[i] ) result_drc = false;
  }
  zAssert( zRawVecDiv (very small or large divisor), result );
  zAssert( zRawVecDivDRC (very small or large divisor), result_drc );
}

/* norm of a vector, the squared norm of which underflows or overflows.
 * a norm smaller than DBL_MIN is a subnormal number, the precision of which is the smallest subnormal
 * number 2^-1074 instead of DBL_EPSILON relative to itself. TEST_NORM_ULP_SUBNORMAL(e) is the smallest
 * subnormal number scaled by 2^-e, which is added to the tolerance of the norm scaled by 2^-e. */
#define N_NORM_TEST 1000
#define TEST_NORM_ULP_SUBNORMAL(e) ldexp( 1.0, -1074-(e) )
void assert_raw_vec_norm_range(void)
{
  const int scale[] = { -1020, -1000, -600, -540, -520, -300, 0, 300, 500, 511, 520, 600, 1000, 1021 };
  double v[TEST_VEC_SIZE], w[TEST_VEC_SIZE], norm, ref;
  const double tiny[] = { 0, 4.9406564584124654e-324, 0 }; /* the smallest subnormal number */
  const double huge[] = { 0, DBL_MAX, 0 };
  const double zero[] = { 0, 0, 0 };
  const double v1[] = { 0.75, -3.5, 2 };
  const double v2[] = { 0.5, -0.25, 0.999 };
  const double inf[] = { 1, HUGE_VAL, 1 };
  double nan[] = { 1, 0, 1 };
  int i, j, k, size;
  bool result = true, result_special;

  for( i=0; i<N_NORM_TEST; i++ ){
    size = zRandI( 1, TEST_VEC_SIZE );
    for( j=0; j<size; j++ ) v[j] = zRandF( -1, 1 );
    for( k=0; k<(int)( sizeof(scale)/sizeof(int) ); k++ ){
      for( j=0; j<size; j++ ) w[j] = ldexp( v[j], scale[k] );
      norm = ldexp( zRawVecNorm( w, size ), -scale[k] );
      for( j=0; j<size; j++ ) w[j] = ldexp( w[j], -scale[k] ); /* exact, even for subnormal components */
      ref = sqrt( zRawVecSqrNorm( w, size ) );
      if( !( fabs( norm - ref ) <= 2 * ( size + 1 ) * DBL_EPSILON * ref + TEST_NORM_ULP_SUBNORMAL(scale[k]) ) ){
        eprintf( "scale 2^%d: norm = %.17g, expected %.17g\n", scale[k], norm, ref );
        result = false;
      }
    }
  }
  nan[1] = sqrt( -zRandF( 1, 2 ) );
  result_special =
    zRawVecNorm( tiny, 3 ) == tiny[1] && zRawVecNorm( huge, 3 ) == DBL_MAX &&
    zRawVecNorm( zero, 3 ) == 0 && zRawVecNorm( zero, 0 ) == 0 &&
    zRawVecNorm( inf, 3 ) == HUGE_VAL && zIsNan( zRawVecNorm( nan, 3 ) ) &&
    zRawVecSqrNorm( huge, 3 ) == HUGE_VAL; /* not a not-a-number by compensated summation */
  zAssert( zRawVecAbsMaxPow2,
    zRawVecAbsMaxPow2( v1, 3 ) == 2 && zRawVecAbsMaxPow2( v2, 3 ) == 0.5 &&
    zRawVecAbsMaxPow2( tiny, 3 ) == tiny[1] && zRawVecAbsMaxPow2( huge, 3 ) == ldexp( 1.0, 1023 ) &&
    zRawVecAbsMaxPow2( zero, 3 ) == 0 && zRawVecAbsMaxPow2( zero, 0 ) == 0 &&
    zRawVecAbsMaxPow2( inf, 3 ) == HUGE_VAL && zIsNan( zRawVecAbsMaxPow2( nan, 3 ) ) );
  zAssert( zRawVecNorm (underflow and overflow), result );
  zAssert( zRawVecNorm (special values), result_special );
}

/* distance and weighted norm of vectors, the squared values of which underflow or overflow. */
void assert_raw_vec_dist_wnorm_range(void)
{
  const int scale[] = { -1020, -600, -520, 0, 520, 600, 1000 };
  const int wscale[] = { -600, 0, 600 }; /* even exponents for exact square roots */
  double v1[TEST_VEC_SIZE], v2[TEST_VEC_SIZE], w[TEST_VEC_SIZE], s1[TEST_VEC_SIZE], s2[TEST_VEC_SIZE], sw[TEST_VEC_SIZE];
  double val, ref;
  const double nan_vec[] = { 1, 0, 1 }, zero[] = { 0, 0, 0 };
  double nan_val;
  int i, j, k, kw, size;
  bool result_dist = true, result_wnorm = true, result_special;

  for( i=0; i<N_NORM_TEST; i++ ){
    size = zRandI( 1, TEST_VEC_SIZE );
    for( j=0; j<size; j++ ){
      v1[j] = zRandF( -1, 1 );
      v2[j] = zRandF( -1, 1 );
      w[j] = zRandF( 0.1, 1 );
    }
    for( k=0; k<(int)( sizeof(scale)/sizeof(int) ); k++ ){
      for( j=0; j<size; j++ ){
        s1[j] = ldexp( v1[j], scale[k] );
        s2[j] = ldexp( v2[j], scale[k] );
      }
      val = ldexp( zRawVecDist( s1, s2, size ), -scale[k] );
      for( j=0; j<size; j++ ){ /* exact, even for subnormal components */
        s1[j] = ldexp( s1[j], -scale[k] );
        s2[j] = ldexp( s2[j], -scale[k] );
      }
      ref = sqrt( zRawVecSqrDist( s1, s2, size ) );
      if( !( fabs( val - ref ) <= 4 * ( size + 1 ) * DBL_EPSILON * ref + TEST_NORM_ULP_SUBNORMAL(scale[k]) ) ){
        eprintf( "zRawVecDist (scale 2^%d): %.17g, expected %.17g\n", scale[k], val, ref );
        result_dist = false;
      }
      for( kw=0; kw<(int)( sizeof(wscale)/sizeof(int) ); kw++ ){
        if( abs( scale[k] + wscale[kw]/2 ) > 1000 ) continue; /* the weighted norm is not representable */
        for( j=0; j<size; j++ ){
          s1[j] = ldexp( v1[j], scale[k] );
          sw[j] = ldexp( w[j], wscale[kw] );
        }
        val = ldexp( zRawVecWNorm( s1, sw, size ), -scale[k]-wscale[kw]/2 );
        for( j=0; j<size; j++ ) s1[j] = ldexp( s1[j], -scale[k] );
        ref = sqrt( zRawVecWSqrNorm( s1, w, size ) );
        if( !( fabs( val - ref ) <= 4 * ( size + 1 ) * DBL_EPSILON * ref + TEST_NORM_ULP_SUBNORMAL(scale[k]+wscale[kw]/2) ) ){
          eprintf( "zRawVecWNorm (scale 2^%d, weight 2^%d): %.17g, expected %.17g\n", scale[k], wscale[kw], val, ref );
          result_wnorm = false;
        }
      }
    }
  }
  nan_val = sqrt( -zRandF( 1, 2 ) );
  memcpy( s1, nan_vec, sizeof(nan_vec) );
  s1[1] = nan_val;
  result_special = zIsNan( zRawVecDist( s1, zero, 3 ) ) && zIsNan( zRawVecWNorm( s1, nan_vec, 3 ) ) &&
    zRawVecDist( zero, zero, 3 ) == 0 && zRawVecWNorm( zero, nan_vec, 3 ) == 0;
  zAssert( zRawVecDist (underflow and overflow), result_dist );
  zAssert( zRawVecWNorm (underflow and overflow), result_wnorm );
  zAssert( zRawVecDist + zRawVecWNorm (special values), result_special );
}

int main(void)
{
  zRandInit();
  assert_raw_vec_get_put();
  assert_raw_vec_equal();
  assert_raw_vec_match();
  assert_raw_vec_misc();
  assert_raw_vec_arith();
  assert_raw_vec_innerprod();
  assert_raw_vec_normalize();
  assert_raw_vec_normalize_tol();
  assert_raw_vec_div_range();
  assert_raw_vec_norm_range();
  assert_raw_vec_dist_wnorm_range();
  return EXIT_SUCCESS;
}
