#include <zm/zm_sf.h>

int main(void)
{
  register int i, j;
  double z, w;

  for( i=1; i<=30; i++ ){
    z = ((double)i) / 10.0;
    for( j=1; j<=50; j++ ){
      w = ((double)j) / 10.0;
      printf( "%.10f %.10f %.10f\n", z, w, zBeta(z,w) );
    }
  }
  return 0;
}
