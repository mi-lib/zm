#include <zm/zm_misc.h>

int main(void)
{
  double x, y;

  for( x=0; x<=30; x+=1.0 ){
    y = zLine( x, 1,-2, 20, 4 );
    printf( "%f %f\n", x, y );
  }
  return 0;
}
