#include <zm/zm_sf.h>

#define N 1000
#define DIV 6

int main(int argc, char *argv[])
{
  int i, k;
  double x, s, c, a;

  for( i=0; i<=N; i++ ){
    x = (double)i/N;
    printf( "%.10g", x );
    for( k=-DIV; k<=DIV; k++ ){
      a = zPI*k/DIV;
      zFresnelIntgGen( x, 0, zPI/3-a, a, &s, &c );
      printf( " %.10g %.10g", c, s );
    }
    printf( "\n" );
  }
  return 0;
}
