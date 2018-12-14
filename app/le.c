#include <zm/zm.h>

void le_usage(char *cmd)
{
  ZRUNERROR( "Usage: %s <eq. file>", cmd );
  exit( 1 );
}

int main(int argc, char *argv[])
{
  zMat a;
  zVec b = NULL, ans = NULL;
  FILE *fp;

  if( argc < 2 ) le_usage( argv[0] );
  if( !( fp = fopen( argv[1], "r" ) ) ){
    ZOPENERROR( argv[1] );
    return 1;
  }
  if( !( a = zMatFRead( fp ) ) ||
      !( b = zVecFRead( fp ) ) ||
      !( ans = zVecAlloc( zMatColSize(a) ) ) ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  zLESolveGauss( a, b, ans );
  zVecWrite( ans );

 TERMINATE:
  zMatFree( a );
  zVecFree( b );
  zVecFree( ans );
  fclose( fp );
  return 0;
}
