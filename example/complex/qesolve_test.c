#include <zm/zm_complex.h>

void heboQEsolve(double a, double b, double c, zComplex ans[])
{
  double d;

  d = b*b - 4*a*c;
  a *= 2;
  if( d >= 0 ){
    d = sqrt( d );
    zComplexCreate( &ans[0], (-b+d)/a, 0 );
    zComplexCreate( &ans[1], (-b-d)/a, 0 );
  } else{
    d = sqrt( -d );
    zComplexCreate( &ans[0], -b/a, d/a );
    zComplexCreate( &ans[1], -b/a,-d/a );
  }
}

void write_ans(zComplex ans[])
{
  printf( " -> x = " );
  zComplexWrite( &ans[0] );
  printf( ", " );
  zComplexWrite( &ans[1] );
  printf( "\n" );
}

void examine(double a, zComplex ans[])
{
  zComplex p;

  printf( "<examine>\n" );
  printf( "%f x^2 + ", a );
  zComplexAdd( &ans[0], &ans[1], &p );
  zComplexMul( &p, -a, &p );
  zComplexWrite( &p );
  printf( " x + " );
  zComplexCMul( &ans[0], &ans[1], &p );
  zComplexMul( &p, a, &p );
  zComplexWrite( &p );
  printf( " = 0\n\n" );
}

void qe_solve(double a, double b, double c)
{
  zComplex ans[2];

  printf( "%f x^2 + %f x + %f = 0\n", a, b, c );

  printf( "<fine>\n" );
  zQESolve( a, b, c, ans );
  write_ans( ans );
  examine( a, ans );

  printf( "<hebo>\n" );
  heboQEsolve( a, b, c, ans );
  write_ans( ans );
  examine( a, ans );
}

int main(void)
{
  qe_solve( 1, 2, 1 );
  qe_solve( 1,-1,-6 );
  qe_solve( 1, 2, 3 );
  qe_solve( 0.000001, 2, 0.000001 );
  return 0;
}
