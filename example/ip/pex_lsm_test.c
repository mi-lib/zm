#include <zm/zm_ip.h>

int main(int argc, char *argv[])
{
  zPexIP pc;
  zVec tvec, xvec;
  double t;

  tvec = zVecCreateList( 6, 0.5, 1.0, 2.0, 3.0, 4.0, 4.5 );
  xvec = zVecCreateList( 6, 1.0, 4.0, 0.0, 3.0, 2.0, 0.0 );
  if( !zPexIPCreateLSM( &pc, 5.0, 4, tvec, xvec ) ){
    ZRUNERROR( "cannot create a polynomial function" );
    zVecFree( tvec );
    zVecFree( xvec );
    return EXIT_FAILURE;
  }

  for( t=0; t<=zPexIPTerm(&pc); t+=0.01 )
    printf( "%.10f %.10f\n", t, zPexIPVal( &pc, t ) );

  zPexIPFree( &pc );
  zVecFree( tvec );
  zVecFree( xvec );
  return 0;
}
