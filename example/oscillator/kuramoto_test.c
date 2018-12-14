#include <zm/zm_oscil.h>

#define DT 0.001
#define T  10.0

void output(double t, zOsc *osc, double x)
{
  printf( "%f %.10f %.10f %.10f %.10f\n", t, zOscOutput(osc), x, sin(zOscOutput(osc)), sin(x) );
}

int main(int argc, char *argv[])
{
  zOsc osc;
  double x;
  register int i, n;

  n = T / DT;
  if( !zOscCreateKura( &osc, zPI, 100, zPI/2 ) ) return 1;
  zOscInit( &osc, 0, 0 );
  output( 0, &osc, 0 );
  for( i=0; i<=n; i++ ){
    x = zPI_2*sin(2*zPI*i*DT);
    zOscUpdate( &osc, x, DT );
    output( i*DT, &osc, x );
  }
  zOscDestroy( &osc );
  return 0;
}
