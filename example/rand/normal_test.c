#include <zm/zm_rand.h>

#define N 80000
#define DIV 1000
int count[DIV];

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  double x;

  zRandInitMT( NULL );
  fp = fopen( "nd", "w" );
  for( i=0; i<N; i++ ){
    x = zRandNormal( NULL );
    fprintf( fp, "%.15f %d\n", x, ++count[(int)((x+5.0)/10.0*DIV)] );
  }
  fclose( fp );
  return 0;
}
