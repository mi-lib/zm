#include <zm/zm_rand.h>

#define N 10000
#define M 100
#define DIV 1000
int count[DIV];

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  double x, n;

  zRandInitMT( NULL );
  n = argc > 1 ? atof(argv[1]) : M;
  fp = fopen( "pd", "w" );
  for( i=0; i<N; i++ ){
    x = zRandPower( NULL, n );
    fprintf( fp, "%.15f %d\n", x, ++count[(int)(x*DIV)] );
  }
  fclose( fp );
  return 0;
}
