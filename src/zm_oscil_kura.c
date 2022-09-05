/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_oscil_kura - nonlinear oscillator: Kuramoto's oscillator.
 */

#include <zm/zm_oscil.h>

typedef struct{
  double wn; /* natural frequency */
  double ke; /* entrainment gain */
  double po; /* phase offset */
  double p;  /* phase */
  double w;  /* state */
} zOscKura;

/* initialize Kuramoto's oscillator. */
static void _zOscInitKura(void *prp, double p0, double w0)
{
  ((zOscKura *)prp)->p = p0;
  ((zOscKura *)prp)->w = w0;
}

/* update the internal state of Kuramoto's oscillator. */
static double _zOscUpdateKura(void *prp, double u, double dt)
{
  zOscKura *kura;

  kura = (zOscKura *)prp;
  kura->w = kura->wn + kura->ke * sin( u - kura->p - kura->po );
  return kura->p += kura->w * dt;
}

/* output a value of Kuramoto's oscillator. */
static double _zOscOutputKura(void *prp)
{
  return ((zOscKura *)prp)->p;
}

static zOscCom _z_osc_com_kura = {
  _zOscInitKura,
  _zOscUpdateKura,
  _zOscOutputKura,
};

/* create Kuramoto's oscillator. */
zOsc *zOscCreateKura(zOsc *osc, double wn, double ke, double po)
{
  zOscKura *kura;

  if( !( kura = zAlloc( zOscKura, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  kura->wn = wn;
  kura->ke = ke;
  kura->po = po;
  osc->prp = kura;
  osc->com = &_z_osc_com_kura;
  return osc;
}
