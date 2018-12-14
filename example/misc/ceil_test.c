#include <zm/zm_misc.h>

int main(void)
{
  register int i;
  double v;

  for( i=-45; i<=45; i++ ){
    v = 0.1 * i;
    printf( "%g %g %g %g %g %g %g\n", v, ceil(v), floor(v), round(v), trunc(v), zCeil(v), zRound(v) );
  }
  return 0;
}
