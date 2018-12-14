#include <zm/zm_rand.h>

#define N 10000
#define DIV 100
int count[DIV];

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  double x, a, b;

  zRandInitMT( NULL );
  a = argc > 1 ? atof( argv[1] ) : 0.5;
  b = argc > 2 ? atof( argv[2] ) : 0.5;
  fp = fopen( "bd", "w" );
  for( i=0; i<N; i++ ){
    x = zRandBeta( NULL, a, b );
    fprintf( fp, "%.15f %d\n", x, ++count[(int)(x*DIV)] );
  }
  fclose( fp );
  return 0;
}
