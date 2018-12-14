#include <zm/zm_pex.h>

int main(void)
{
  double x;
  zPex p;

  p = zPexCreateList( 3, 1.0,-1.0, 3.0,-1.0 );
  for( x=-2.0; x<=3.0; x+=0.01 )
    printf( "%.10g\n", zPexVal( p, x ) );
  zPexFree( p );
  return 0;
}
