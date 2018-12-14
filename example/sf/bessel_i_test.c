#include <zm/zm_sf.h>

int main(void)
{
  register int i;
  double x;

  for( i=-1000; i<=1000; i++ ){
    x = 0.02 * ((double)i);
    printf( "%.16f %.16f %.16f %.16f %.16f\n", x, zBesselI(0,x), zBesselI(1,x), zBesselI(2,x), zBesselI(3,x) );
  }
  return 0;
}
