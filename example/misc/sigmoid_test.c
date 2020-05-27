#include <zm/zm_misc.h>

int main(void)
{
  double x, y;

  for( x=-30; x<=30; x+=1 ){
    y = 5 * zSigmoid( 0.25*x );
    printf( "%f %f\n", x, y );
  }
  return 0;
}
