/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_intg - numerical integration.
 */

#include <zm/zm_intg.h>

/* zIntgInit
 * - initialize integrator.
 */
void zIntgInit(zIntg *intg, double s0, double x0)
{
  intg->s = s0;
  intg->x0 = x0;

  intg->_dt = 0;
  intg->_x = 0;
}

/* zIntgRect
 * - rectangular integration.
 */
double zIntgRect(zIntg *intg, double x, double dt)
{
  intg->s += intg->x0 * dt;
  intg->x0 = x;
  return intg->s;
}

/* zIntgTR
 * - trapezoidal integration.
 */
double zIntgTR(zIntg *intg, double x, double dt)
{
  intg->s += 0.5 * ( intg->x0 + x ) * dt;
  intg->x0 = x;
  return intg->s;
}

/* zIntgQuad
 * - integration with an arrangement of Simpson's formula.
 */
double zIntgQuad(zIntg *intg, double x, double dt)
{
  double a, b, c, d;

  if( intg->_dt == 0 ){
    intg->_dt = dt;
    intg->_x = intg->x0;
    return zIntgTR( intg, x, dt );
  }
  d = dt + intg->_dt;
  a = dt / ( 6*intg->_dt ) + 0.5;
  b = ( dt/3 + 0.5*intg->_dt ) / d;
  c = - dt*dt / ( 6*intg->_dt*d );
  intg->s += ( intg->x0 * a + x * b + intg->_x * c ) * dt;
  intg->_x = intg->x0;
  intg->_dt = dt;
  intg->x0 = x;
  return intg->s;
}
