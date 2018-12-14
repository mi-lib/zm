#include <zm/zm_rand.h>

int main(void)
{
  double x, y, z;
  long long n;
  FILE *fp1, *fp2;
  register int i;

  n = 2 << 10;
  printf( "%lld\n", n );
  zRandInit();
  zRandInitMT( NULL );
  fp1 = fopen( "nm", "w" );
  fp2 = fopen( "mt", "w" );
  for( i=0; i<n; i++ ){
    x = zRandF( 0, 1 );
    y = zRandF( 0, 1 );
    z = zRandF( 0, 1 );
    fprintf( fp1, "%.15f %.15f %.15f\n", x, y, z );
    x = zRandMTF( NULL, 0, 1 );
    y = zRandMTF( NULL, 0, 1 );
    z = zRandMTF( NULL, 0, 1 );
    fprintf( fp2, "%.15f %.15f %.15f\n", x, y, z );
  }
  fclose( fp1 );
  fclose( fp2 );
  return 0;
}
