#include <zm/zm_rand.h>

#define N 1000

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  double x, min, max;
  int d, div, *count;

  zRandInitMT( NULL );
  min = argc > 1 ? atof( argv[1] ) : -100000;
  max = argc > 2 ? atof( argv[2] ) :  100000;
  div = max - min + 1;
  count = calloc( div, sizeof(int) );
  fp = fopen( "cd", "w" );
  for( i=0; i<N; i++ ){
    x = zRandCauchy( NULL );
    d = (int)x > max ? div-1 : ( (int)x < min ? 0 : (int)x + div/2 );
    fprintf( fp, "%.15f %d\n", x, ++count[d] );
  }
  fclose( fp );
  return 0;
}
