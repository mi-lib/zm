#include <zm/zm_rand.h>

#define N 10000
#define M 1000
#define DIV M
int count[DIV];

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  int x;
  double p;

  zRandInitMT( NULL );
  p = argc > 1 ? atof(argv[1]) : 0.5;
  fp = fopen( "bd", "w" );
  for( i=0; i<N; i++ ){
    x = zRandBinom( NULL, M, p );
    fprintf( fp, "%d %d\n", x, ++count[x] );
  }
  fclose( fp );
  return 0;
}
