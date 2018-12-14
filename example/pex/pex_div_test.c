#include <zm/zm_pex.h>

int main(void)
{
  zPex p, f, q, r, g;

  p = zPexCreateList( 5, 1.0,-3.0, 7.0, 3.0,-4.0, 2.0 );
  f = zPexCreateList( 3, 3.0, 2.0, 2.0, 2.0 );

  zPexDiv( p, f, &q, &r );
  g = zPexMul( f, q );
  zPexAddDRC( g, r );

  printf( "P(x) = " );
  zPexExprX( p );
  printf( "f(x) = " );
  zPexExprX( f );
  printf( "Q(x) = " );
  zPexExprX( q );
  printf( "R(x) = " );
  zPexExprX( r );

  zEndl();
  printf( "f(x)Q(x) + R(x) = " );
  zPexExprX( g );
  printf( "%s\n", zBoolExpr( zPexIsEqual( p, g ) ) );

  zPexFree( p );
  zPexFree( f );
  zPexFree( q );
  zPexFree( r );
  return 0;
}
