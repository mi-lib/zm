#include <zm/zm_pex.h>

int main(void)
{
  zPex p, pd, pi;

  p = zPexCreateList( 5, 1.0,-3.0, 7.0, 3.0,-4.0, 2.0 );

  pd = zPexDif( p );
  pi = zPexIntg( p );

  printf( "f(x) = " );
  zPexExprX( p );
  printf( "f'(x) = " );
  zPexExprX( pd );
  printf( "F(x) = " );
  zPexExprX( pi );

  zPexFree( p );
  zPexFree( pd );
  zPexFree( pi );
  return 0;
}
