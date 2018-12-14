#include <zm/zm_pex.h>

int main(void)
{
  zPex p;
  FILE *fp;

  if( !( fp = fopen( "pex_test.zpx", "r" ) ) )
    return 1;
  p = zPexFRead( fp );
  fclose( fp );

  zPexExprX( p );
  zPexWrite( p );
  zPexFree( p );
  return 0;
}
