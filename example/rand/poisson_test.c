#include <zm/zm_rand.h>

#define N 10000
#define DIV 1000
int count[DIV];

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  int x;
  double alpha;

  zRandInitMT( NULL );
  alpha = argc > 1 ? atof( argv[1] ) : 10;
  fp = fopen( "pd", "w" );
  for( i=0; i<N; i++ ){
    x = zRandPoisson( NULL, alpha );
    if( x >= DIV ) x = DIV - 1;
    fprintf( fp, "%d %d\n", x, ++count[x] );
  }
  fclose( fp );
  return 0;
}
