/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_oscil_vdp - nonlinear oscillator: van der Pol's oscillator.
 */

#include <zm/zm_oscil.h>

typedef struct{
  double omega;  /* natural angular frequency */
  double eps;    /* damping coefficient */
  double amp;    /* output amplitude */
  double _x, _v; /* internal state */
} zOscVDP;

/* initialize van der Pol's oscillator. */
static void _zOscInitVDP(void *prp, double x0, double v0)
{
  ((zOscVDP *)prp)->_x = x0 / ((zOscVDP *)prp)->amp;
  ((zOscVDP *)prp)->_v = v0 / ((zOscVDP *)prp)->amp;
}

/* output a value of van der Pol's oscillator. */
static double _zOscOutputVDP(void *prp)
{
  return ((zOscVDP *)prp)->amp * ((zOscVDP *)prp)->_x;
}

/* update the internal state of van der Pol's oscillator. */
static double _zOscUpdateVDP(void *prp, double u, double dt)
{
  zOscVDP *vdp;
  double a;

  vdp = (zOscVDP *)prp;
  a = -vdp->omega * ( vdp->eps*( zSqr(vdp->_x) - 1 ) * vdp->_v + vdp->omega*vdp->_x ) - u;
  vdp->_v += a       * dt;
  vdp->_x += vdp->_v * dt;
  return _zOscOutputVDP( prp );
}

static zOscCom _z_osc_com_vdp = {
  _zOscInitVDP,
  _zOscUpdateVDP,
  _zOscOutputVDP,
};

/* create van der Pol's oscillator. */
zOsc *zOscCreateVDP(zOsc *osc, double t, double damp, double amp)
{
  zOscVDP *vdp;

  if( zIsTiny(t) ){
    ZRUNERROR( ZM_ERR_OSCIL_INVALID_TERM );
    return NULL;
  }
  if( !( vdp = zAlloc( zOscVDP, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  vdp->omega = zPIx2 / t;
  vdp->eps   = damp;
  vdp->amp   = amp * 0.5; /* natural amplitude of VDP is 2 */
  osc->prp = vdp;
  osc->com = &_z_osc_com_vdp;
  return osc;
}
