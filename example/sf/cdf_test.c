#include <zm/zm_sf.h>

#define W 5.0
#define DIV 1000

int main(int argc, char *argv[])
{
  int i;
  double x, ef, cf;

  for( i=0; i<=DIV; i++ ){
    x = W * ( (double)i/DIV - 0.5 );
    ef = zErf(x);
    cf = zCDF(x);
    printf( "%.10g %.10g %.10g\n", x, ef, cf );
  }
  return 0;
}
