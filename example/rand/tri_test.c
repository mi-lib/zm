#include <zm/zm_rand.h>

#define N 100000
#define DIV 1000
int count[2*DIV+1];

int main(int argc, char *argv[])
{
  FILE *fp;
  register int i;
  double x;

  zRandInitMT( NULL );
  fp = fopen( "td", "w" );
  for( i=0; i<N; i++ ){
    x = zRandTri( NULL );
    fprintf( fp, "%.15f %d\n", x, ++count[(int)((x+1)*DIV)] );
  }
  fclose( fp );
  return 0;
}
