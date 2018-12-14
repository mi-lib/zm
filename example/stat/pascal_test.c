#include <zm/zm_stat.h>

double _zCombi(int n, int i)
{
  if( n <= 0 ) return 1;
  if( i == 0 || i == n ) return 1;
  if( i < 0 || i > n ) return NAN;
  return _zCombi( n-1, i-1 ) + _zCombi( n-1, i );
}

#define MAX 100

int main(int argc, char *argv[])
{
  int n, i, j;
  double c[MAX];

  n = argc > 1 ? atoi( argv[1] ) : 6;
  printf( "[Pascal's triangle]\n" );
  printf( "*** naive method ***\n" );
  for( i=0; i<=n; i++ ){
    for( j=0; j<=i; j++ )
      printf( " %d", (int)zCombi(i,j) );
    printf( "\n" );
  }
  printf( "*** recursive method ***\n" );
  for( i=0; i<=n; i++ ){
    for( j=0; j<=i; j++ )
      printf( " %d", (int)_zCombi(i,j) );
    printf( "\n" );
  }
  printf( "*** batch method ***\n" );
  for( i=0; i<=n; i++ ){
    zCombiSeries( i, MAX, c );
    for( j=0; j<=i; j++ )
      printf( " %d", (int)c[j] );
    printf( "\n" );
  }
  return 0;
}
