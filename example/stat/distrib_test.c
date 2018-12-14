/* distribution test */
/* this example shows the use of
 *  1: normal distribution
 *  2: poisson's distribution
 *  3: binomial distribution
 */
#include <zm/zm_stat.h>

#define DX (double)1.0

int main(void)
{
  register int i;
  double x;

  for( i=0; i<=20; i++ ){
    x = DX * i;
    printf( "%f %f %f %f\n", x, zNormalDistrib( x, 3, 2 ), zPoissonDistrib( x, 10 ), zBinDistrib( x, 20, 0.5 ) );
  }
  return 0;
}
