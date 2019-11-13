#include <zm/zm_ip.h>

int main(int argc, char *argv[])
{
  zPexIP pc;
  zVec tvec, xvec;
  double t;

  tvec = zVecCreateList( 3, 1.0, 2.0, 4.0 );
  xvec = zVecCreateList( 3, 4.0, 0.0, 2.0 );
  zPexIPCreateLSM( &pc, 5.0, 4, tvec, xvec );

  for( t=0; t<=zPexIPTerm(&pc); t+=0.01 )
    printf( "%.10f %.10f\n", t, zPexIPVal( &pc, t ) );

  zPexIPFree( &pc );
  return 0;
}
