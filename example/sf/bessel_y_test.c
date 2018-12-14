#include <zm/zm_sf.h>

int main(void)
{
  register int i;
  double x;

  for( i=1; i<=1000; i++ ){
    x = 0.02 * ((double)i);
    printf( "%.16f %.16f %.16f %.16f %.16f\n", x, zBesselY(0,x), zBesselY(1,x), zBesselY(2,x), zBesselY(3,x) );
  }
  return 0;
}
