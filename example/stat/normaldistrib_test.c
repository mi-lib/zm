#include <zm/zm_stat.h>

#define DX (double)1.0

int main(int argc, char *argv[])
{
  register int i;
  double x;

  for( i=-17; i<=20; i++ ){
    x = DX * i;
    printf( "%.10g %.10g %.10g\n", x, zNormalDistrib( x, 3, 2 ), zNormalCumDistrib( x, 3, 2 ) );
  }
  return 0;
}
