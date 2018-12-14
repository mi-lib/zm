#include <zm/zm_rand.h>

#define N 10000
#define DIV 1000
int count[DIV];

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  double x;

  zRandInitMT( NULL );
  fp = fopen( "ed", "w" );
  for( i=0; i<N; i++ ){
    x = zRandExp( NULL );
    fprintf( fp, "%.15f %d\n", x, ++count[(int)(x*0.1*DIV)] );
  }
  fclose( fp );
  return 0;
}
