#include <zm/zm_misc.h>

#define DIV 100

int main(void)
{
  double x, y, t;
  register int i;

  for( i=0; i<=DIV; i++ ){
    t = (double)i / DIV;
    x = 7 * zCycloidX( t ) + 5;
    y = 3 * zCycloidY( t ) + 4;
    printf( "%f %f\n", x, y );
  }
  return 0;
}
