#include <zm/zm_rand.h>

#define N 10000
#define DIV 1000
int count[DIV];

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  double x;
  double alpha;

  zRandInitMT( NULL );
  alpha = argc > 1 ? atof(argv[1]) : 0.5;
  fp = fopen( "wd", "w" );
  for( i=0; i<N; i++ ){
    x = zRandWeibull( NULL, alpha );
    fprintf( fp, "%.15f %d\n", x, ++count[(int)x>=DIV ? DIV-1 : (int)x] );
  }
  fclose( fp );
  return 0;
}
