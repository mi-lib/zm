#include <zm/zm_rand.h>

#define N 10000
int *count;

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  int div, x;
  double p;

  zRandInitMT( NULL );
  p = argc > 1 ? atof(argv[1]) : 0.5;
  div = 10.0 / p;
  count = calloc( div, sizeof(int) );
  fp = fopen( "gd", "w" );
  for( i=0; i<N; i++ ){
    x = zRandGeo( NULL, p );
    fprintf( fp, "%d %d\n", x, ++count[x] );
  }
  fclose( fp );
  return 0;
}
