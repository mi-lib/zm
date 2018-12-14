#include <zm/zm_rand.h>

#define N 10000
#define DIV 1000
int count[DIV];

int main(int argc, char *argv[])
{
  double x;
  FILE *fp;
  register int i;
  double n;
  int d;

  zRandInitMT( NULL );
  n = argc > 1 ? atof(argv[1]) : 2.0;
  if( n < 2.0 ){
    ZRUNWARN( "probably unappropriate n=%f, modified to be 2.0", n );
    n = 2.0;
  }
  fp = fopen( "td", "w" );
  for( i=0; i<N; i++ ){
    x = zRandT( NULL, n );
    d = (x+10)/20 * DIV;
    fprintf( fp, "%.15f %d\n", x, ++count[d] );
  }
  fclose( fp );
  return 0;
}
