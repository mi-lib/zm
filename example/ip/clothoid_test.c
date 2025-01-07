#include <zm/zm_ip.h>

#define N 100

void plot(zClothoid *cl)
{
  int i;
  double s, x, y;

  for( i=0; i<=N; i++ ){
    zClothoidXY( cl, ( s = (double)i/N ), &x, &y );
    printf( "%.10g %.10g %.10g\n", s, x, y );
  }
}

int main(int argc, char *argv[])
{
  zClothoid cl;

  zClothoidCreateSegment( &cl, 0, 0,-0.5*zPI, 2, 1, 0.5*zPI );
  plot( &cl );
  return 0;  
}
