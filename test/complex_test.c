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

int main(int argc, char *argv[])
{
  zRandInit();
  zAssert( zComplexValIsIncluded, assert_isincluded() );
  zAssert( zQESolve, assert_qesolve() );
  zAssert( zCESolve, assert_cesolve() );
  return 0;
}
