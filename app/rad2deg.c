#include <zm/zm.h>

void usage(void)
{
  fprintf( stderr, "Usage: rad2deg <val1> <val2> ...\n" );
  fprintf( stderr, "       rad2deg <file> <id>\n" );
  exit( 1 );
}

void rad2deg_file(FILE *fp, int id)
{
  char buf[BUFSIZ];
  double val;
  int i;

  while( fgets( buf, BUFSIZ, fp ) ){
    val = 0;
    for( i=0; i<id; i++ )
      val = zSDouble( buf );
    printf( "%.10f\n", zRad2Deg(val) );
  }
  fclose( fp );
}

int main(int argc, char *argv[])
{
  FILE *fp;

  if( argc < 2 ) usage();
  if( ( fp = fopen( argv[1], "r" ) ) ){
    if( argc < 3 ) usage();
    rad2deg_file( fp, atoi( argv[2] ) );
  } else{
    while( *++argv )
      printf( "%.10f ", zRad2Deg(atof(*argv)) );
    printf( "\n" );
  }
  return 0;
}
