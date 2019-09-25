#include <zm/zm_pex.h>

int main(void)
{
  zPex p;
  FILE *fp;

  if( !( fp = fopen( "pex_test.dat", "r" ) ) )
    return 1;
  p = zPexFScan( fp );
  fclose( fp );

  zPexExprX( p );
  zPexPrint( p );
  zPexFree( p );
  return 0;
}
