#include <zm/zm_sf.h>

int main(void)
{
  register int i;
  double x;

  for( i=-1000; i<=1000; i++ ){
    x = 0.02 * ((double)i);
    printf( "%.10f %.10f %.10f %.10f %.10f\n", x, zBesselJ(0,x), zBesselJ(1,x), zBesselJ(2,x), zBesselJ(3,x) );
  }
  return 0;
}
