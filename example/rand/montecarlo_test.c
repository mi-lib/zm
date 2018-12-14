#include <zm/zm_rand.h>

#define N 100000

int test_glibc(int n)
{
  double x, y;
  int nc;
  register int i;

  zRandInit();
  for( nc=0, i=0; i<n; i++ ){
    x = zRandF( -1, 1 );
    y = zRandF( -1, 1 );
    if( x*x + y*y <= 1 ) nc++;
  }
  return nc;
}

int test_mt(int n)
{
  double x, y;
  int nc;
  register int i;

  zRandInitMT( NULL );
  for( nc=0, i=0; i<n; i++ ){
    x = zRandMTF( NULL, -1, 1 );
    y = zRandMTF( NULL, -1, 1 );
    if( x*x + y*y <= 1 ) nc++;
  }
  return nc;
}

void eval(int n, int (*test)(int), const char *label)
{
  int nc;
  double pi_e;

  nc = test( n );
  pi_e = 4.0 * nc / n;
  printf( "%s: (estim) %.16g, (error) %.16g\n", label, pi_e, (pi_e-zPI)/zPI );
}

int main(int argc, char *argv[])
{
  int n;

  n = argc > 1 ? atoi( argv[1] ) : N;
  eval( n, test_glibc, "Glibc" );
  eval( n, test_mt,    "MT   " );
  return 0;
}
