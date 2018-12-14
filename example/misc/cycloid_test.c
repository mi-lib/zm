#include <zm/zm_misc.h>

#define DIV 100

int main(void)
{
  double x, y, t;
  register int i;

  for( i=0; i<=DIV; i++ ){
    t = (double)i / DIV;
    x = zCycloidX( 5, 12, t );
    y = zCycloidY( 3,  4, t );
    printf( "%f %f\n", x, y );
  }
  return 0;
}
