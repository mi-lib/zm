#include <zm/zm_misc.h>

#define N 100

int main(void)
{
  double x;
  int i, j;

  for( i=-N; i<=2*N; i++ ){
    x = (double)i / N;
    for( j=0; j<5; j++ ){
      printf( " %g %g", zSmoothStep( x, j ), zSmoothStepDif( x, j ) );
    }
    printf( "\n" );
  }
  return 0;
}
