#include <zm/zm_rand.h>

#define N 10000
#define DIV 30
int count[DIV];

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  double x;

  zRandInitMT( NULL );
  fp = fopen( "ld", "w" );
  for( i=0; i<N; i++ ){
    x = zRandLog( NULL );
    fprintf( fp, "%.15f %d\n", x, ++count[(int)(x+DIV/2)] );
  }
  fclose( fp );
  return 0;
}
