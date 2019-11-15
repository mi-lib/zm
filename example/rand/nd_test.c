#include <zm/zm_rand.h>

#define N 5

int main(void)
{
  long long n;
  FILE *fp;
  char filename[BUFSIZ];
  register int i, j;

  n = 2 << 10;
  zRandInit();
  for( i=0; i<N; i++ ){
    sprintf( filename, "%d", i );
    fp = fopen( filename, "w" );
    for( j=0; j<n; j++ )
      fprintf( fp, "%.15f\n", zRandND(NULL,1.0,0.5*pow(2,i)) );
    fclose( fp );
  }
  return 0;
}
