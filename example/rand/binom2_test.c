#include <zm/zm_rand.h>

#define N 10000
#define DIV 20
int count_x[DIV], count_y[DIV];

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  double x, y;
  double r;

  zRandInitMT( NULL );
  r = argc > 1 ? atof(argv[1]) : 0.5;
  fp = fopen( "2d", "w" );
  for( i=0; i<N; i++ ){
    zRandBinom2( NULL, r, &x, &y );
    fprintf( fp, "%.15f %d %.15f %d\n", x, ++count_x[(int)(x+(double)DIV/2)], y, ++count_y[(int)(y+(double)DIV/2)] );
  }
  fclose( fp );
  return 0;
}
