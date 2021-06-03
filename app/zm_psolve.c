#include <zm/zm.h>

void psolve_usage(char *cmd)
{
  ZRUNERROR( "Usage: %s <pex file>", cmd );
  exit( 1 );
}

int main(int argc, char *argv[])
{
  zPex pex;
  zCVec ans;
  int dim, i;
  FILE *fp;

  if( argc < 2 ){
    fp = stdin;
  } else{
    if( strcmp( argv[1], "-help" ) == 0 )
      psolve_usage( argv[0] );
    if( !( fp = fopen( argv[1], "r" ) ) ){
      ZOPENERROR( argv[1] );
      return EXIT_FAILURE;
    }
  }
  if( !( pex = zPexFScan( fp ) ) ){
    ZALLOCERROR();
    return EXIT_FAILURE;
  }
  if( fp != stdin ) fclose( fp );

  dim = zPexDim( pex );
  if( !( ans = zCVecAlloc( dim ) ) ) return EXIT_FAILURE;
  printf( "Equation:\n" );
  zPexExprX( pex );
  printf( "= 0\n" );
  zPexDKA( pex, ans, zTOL, 0 );
  printf( "Answer(s):\n" );
  for( i=0; i<dim; i++ ){
    zComplexPrint( zCVecElemNC(ans,i) );
    zEndl();
  }
  zCVecFree( ans );
  zPexFree( pex );
  return 0;
}
