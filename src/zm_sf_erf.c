/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_sf_erf - special functions : Gauss's error function.
 */

#include <zm/zm_sf.h>

/* zErfN
 * - error function approximated by Taylor series.
 */
double zErfN(double x, int n)
{
  double x2, s, d;
  register int i, j;

  x2 = x*x;
  for( s=0, i=0; i<=n; i++ ){
    d = 1;
    for( j=1; j<=i; j++ )
      d *= -x2/j;
    s += d*x/(2*i+1);
  }
  return 2*s / sqrt(zPI);
}
