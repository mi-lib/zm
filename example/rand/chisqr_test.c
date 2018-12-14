#include <zm/zm_rand.h>

#define N 10000
#define DIV 1000
int count[DIV];

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  double x, a;

  zRandInitMT( NULL );
  a = argc > 1 ? atof( argv[1] ) : 0.5;
  fp = fopen( "cd", "w" );
  for( i=0; i<N; i++ ){
    x = zRandChiSqr( NULL, a );
    fprintf( fp, "%.15f %d\n", x, ++count[(int)(0.05*x/a*DIV)] );
  }
  fclose( fp );
  return 0;
}
