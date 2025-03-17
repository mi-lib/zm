#include <zm/zm_rand.h>

int main(void)
{
  double x, y, z;
  long long n;
  FILE *fp;
  int i;

  n = 2 << 10;
  printf( "%lld\n", n );
  zRandInit();
  fp = fopen( "mt", "w" );
  for( i=0; i<n; i++ ){
    x = zRandF( 0, 1 );
    y = zRandF( 0, 1 );
    z = zRandF( 0, 1 );
    fprintf( fp, "%.15f %.15f %.15f\n", x, y, z );
  }
  fclose( fp );
  return 0;
}
