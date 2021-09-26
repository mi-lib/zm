#include <zm/zm_sf.h>

#define N 1000

int main(int argc, char *argv[])
{
  int i;
  double x, s1, c1, s2, c2;

  for( i=-N; i<=N; i++ ){
    x = 10*(double)i/N;
    zFresnelIntg( x, &s1, &c1 );
    zFresnelIntgPI_2( x, &s2, &c2 );
    printf( "%.10g %.10g %.10g %.10g %.10g\n", x, s1, c1, s2, c2 );
  }
  return 0;
}
