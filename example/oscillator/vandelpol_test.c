#include <zm/zm_oscil.h>

#define DT 0.005
#define T  10.0

int main(int argc, char *argv[])
{
  zOsc osc;
  double term, eps, amp, x0, v0;
  int i, n;

  n = T / DT;
  term = argc > 1 ? atof(argv[1]) : 1.0;
  eps  = argc > 2 ? atof(argv[2]) : 1.5;
  amp  = argc > 3 ? atof(argv[3]) : 1.0;
  x0   = argc > 4 ? atof(argv[4]) : 0.5;
  v0   = argc > 5 ? atof(argv[5]) : 0.0;
  if( !zOscCreateVDP( &osc, term, eps, amp ) ) return 1;
  zOscInit( &osc, x0, v0 );
  for( i=0; i<=n; i++ ){
    zOscUpdate( &osc, 0, DT );
    printf( "%f %.10f\n", i*DT, zOscOutput(&osc) );
  }
  zOscDestroy( &osc );
  return 0;
}
