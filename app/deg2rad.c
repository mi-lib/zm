#include <zm/zm.h>

void usage(void)
{
  fprintf( stderr, "Usage: deg2rad <val1> <val2> ...\n" );
  fprintf( stderr, "       deg2rad <file> <id>\n" );
  exit( 1 );
}

void deg2rad_file(FILE *fp, int id)
{
  char buf[BUFSIZ];
  double val;
  int i;

  while( fgets( buf, BUFSIZ, fp ) ){
    val = 0;
    for( i=0; i<id; i++ )
      val = zSDouble( buf );
    printf( "%.10f\n", zDeg2Rad(val) );
  }
  fclose( fp );
}

int main(int argc, char *argv[])
{
  FILE *fp;

  if( argc < 2 ) usage();
  if( ( fp = fopen( argv[1], "r" ) ) ){
    if( argc < 3 ) usage();
    deg2rad_file( fp, atoi( argv[2] ) );
  } else{
    while( *++argv )
      printf( "%.10f ", zDeg2Rad(atof(*argv)) );
    printf( "\n" );
  }
  return 0;
}
