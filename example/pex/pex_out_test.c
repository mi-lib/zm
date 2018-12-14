#include <zm/zm_pex.h>

int main(void)
{
  double x;
  zPex p;

  p = zPexCreateList( 3, 1.0,-1.0, 3.0,-1.0 );
  zPexExprX( p );

  for( x=-2.0; x<=3.0; x+=0.01 )
    printf( "%f %f %f %f %f %f\n",
      zPexDifVal( p, 0, x ),
      zPexDifVal( p, 1, x ),
      zPexDifVal( p, 2, x ),
      zPexDifVal( p, 3, x ),
      zPexDifVal( p, 4, x ),
      zPexDifVal( p, 5, x ) );

  zPexFree( p );
  return 0;
}
