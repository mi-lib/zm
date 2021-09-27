#include <zm/zm_sf.h>

#define N 1000
#define DIV 3

int main(int argc, char *argv[])
{
  int i, k;
  double x, s, c;

  for( i=-N; i<=N; i++ ){
    x = 10*(double)i/N;
    zFresnelIntg( x, &s, &c );
    printf( "%.10g %.10g %.10g", x, s, c );
    for( k=1; k<=DIV; k++ ){
      zFresnelIntgScale( x, 1+(zPI_2-1)*(double)k/DIV, &s, &c );
      printf( " %.10g %.10g", s, c );
    }
    zFresnelIntgPI_2( x, &s, &c );
    printf( " %.10g %.10g\n", s, c );
  }
  return 0;
}
