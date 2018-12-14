#include <zm/zm.h>

void pexp_usage(char *cmd)
{
  ZRUNERROR( "Usage: %s <.zv file>\n", cmd );
  ZRUNERROR( "     : %s\n", cmd );
  ZRUNERROR( " (format) <number of factors> <factor1> <factor2> ...\n" );
  ZRUNERROR( "     : %s -help\n", cmd );
  ZRUNERROR( " (to show this message.)\n" );
  exit( EXIT_SUCCESS );
}

int main(int argc, char *argv[])
{
  zPex pex;
  zVec fact;
  FILE *fp;

  if( argc < 2 ){
    fp = stdin;
  } else{
    if( strcmp( argv[1], "-help" ) == 0 )
      pexp_usage( argv[0] );
    if( !( fp = fopen( argv[1], "r" ) ) ){
      ZOPENERROR( argv[1] );
      return EXIT_FAILURE;
    }
  }
  if( !( fact = zVecFRead( fp ) ) ){
    ZALLOCERROR();
    return EXIT_FAILURE;
  }
  if( fp != stdin ) fclose( fp );
  if( !( pex = zPexExp( fact ) ) ){
    ZALLOCERROR();
    return EXIT_FAILURE;
  }

  printf( "Factor(s):\n" );
  zVecWrite( fact );
  printf( "Expr:\n" );
  zPexExprX( pex );
  zVecFree( fact );
  zPexFree( pex );
  return EXIT_SUCCESS;
}
