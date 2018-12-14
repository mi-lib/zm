#include <zm/zm_rand.h>

#define N 10000
#define DIV 10000
int count[DIV];

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  double x, a, b;
  int d;

  zRandInitMT( NULL );
  a = argc > 1 ? atof( argv[1] ) : 0.1;
  b = argc > 2 ? atof( argv[2] ) : 1.0;
  fp = fopen( "fd", "w" );
  for( i=0; i<N; i++ ){
    x = zRandFD( NULL, a, b );
    d = x * DIV;
    if( d >= DIV ) d = DIV-1;
    fprintf( fp, "%.15f %d\n", x, ++count[d] );
  }
  fclose( fp );
  return 0;
}
