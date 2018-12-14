#include <zm/zm_stat.h>

double _zCombi(int n, int i)
{
  if( n <= 0 ) return 1;
  if( i == 0 || i == n ) return 1;
  if( i < 0 || i > n ) return NAN;
  return _zCombi( n-1, i-1 ) + _zCombi( n-1, i );
}

void measure_combi(int n, int i)
{
  clock_t t1, t2;
  double c;

  printf( "*** naive method ***\n" );
  t1 = clock();
  c = zCombi( n, i );
  t2 = clock();
  printf( "%g: +%ld\n", c, (long)(t2-t1) );

  printf( "*** recursive method ***\n" );
  t1 = clock();
  c = _zCombi( n, i );
  t2 = clock();
  printf( "%g: +%ld\n", c, (long)(t2-t1) );
}

#define MAX 5000
void measure_series(int n)
{
  register int i;
  clock_t t1, t2;
  double c[MAX];

  printf( "*** naive method ***\n" );
  t1 = clock();
  for( i=0; i<=n; i++ )
    c[i] = zCombi( n, i );
  t2 = clock();
  for( i=0; i<=n; i++ )
    printf( " %g", c[i] );
  printf( "\n +%ld\n", (long)(t2-t1) );

  printf( "*** recursive method ***\n" );
  t1 = clock();
  zCombiSeries( n, MAX, c );
  t2 = clock();
  for( i=0; i<=n; i++ )
    printf( " %g", c[i] );
  printf( "\n +%ld\n", (long)(t2-t1) );
}

int main(int argc, char *argv[])
{
  int n, i;

  n = argc > 1 ? atoi( argv[1] ) : 20;
  i = argc > 2 ? atoi( argv[2] ) : 1;
  /*
  measure_combi( n, i );
  */
  measure_series( n );
  return 0;
}
