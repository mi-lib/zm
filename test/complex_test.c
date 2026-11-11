#include <zm/zm.h>

#define NUM 100

#define TOL ( 1.0e-7 )

bool assert_qesolve_check(double a, double b, double c, zComplex *ans)
{
  zComplex r1, r2;

  zComplexMul( ans, a, &r1 );
  r1.re += b;
  zComplexCMul( ans, &r1, &r2 );
  r2.re += c;
  return zComplexIsTol( &r2, TOL );
}

bool assert_qesolve(void)
{
  int i;
  double a, b, c;
  zComplex ans[2];

  for( i=0; i<NUM; i++ ){
    a = zRandF(-3,3); if( a == 0 ) a = 1;
    b = zRandF(-3,3);
    c = zRandF(-3,3);
    zQESolve( a, b, c, ans );
    if( !assert_qesolve_check( a, b, c, &ans[0] ) ||
        !assert_qesolve_check( a, b, c, &ans[1] ) ) return false;
  }
  return true;
}

bool assert_cesolve_check(double a, double b, double c, double d, zComplex *ans)
{
  zComplex r1, r2, r3;

  zComplexMul( ans, a, &r1 );
  r1.re += b;
  zComplexCMul( ans, &r1, &r2 );
  r2.re += c;
  zComplexCMul( ans, &r2, &r3 );
  r3.re += d;
  return zComplexIsTol( &r3, TOL );
}

bool assert_cesolve(void)
{
  int i;
  double a, b, c, d;
  zComplex ans[3];

  for( i=0; i<NUM; i++ ){
    a = zRandF(-3,3); if( a == 0 ) a = 1;
    b = zRandF(-3,3);
    c = zRandF(-3,3);
    d = zRandF(-3,3);
    zCESolve( a, b, c, d, ans );
    if( !assert_cesolve_check( a, b, c, d, &ans[0] ) ||
        !assert_cesolve_check( a, b, c, d, &ans[1] ) ||
        !assert_cesolve_check( a, b, c, d, &ans[2] ) ) return false;
  }
  return true;
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

int main(int argc, char *argv[])
{
  zRandInit();
  zAssert( zComplexDivDRC + _zComplexDivDRC, assert_div_drc() );
  zAssert( zComplexCMul + zComplexCMulConj (in-place), assert_cmul_inplace() );
  zAssert( zComplexValIsIncluded, assert_isincluded() );
  zAssert( zQESolve, assert_qesolve() );
  zAssert( zCESolve, assert_cesolve() );
  return 0;
}
