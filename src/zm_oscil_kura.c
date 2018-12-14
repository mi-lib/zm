/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_oscil_kura - nonlinear oscillator: Kuramoto's oscillator.
 */

#include <zm/zm_oscil.h>

static void _zOscInitKura(void *prp, double p0, double w0);
static double _zOscUpdateKura(void *prp, double u, double dt);
static double _zOscOutputKura(void *prp);

typedef struct{
  double wn; /* natural frequency */
  double ke; /* entrainment gain */
  double po; /* phase offset */
  double p;  /* phase */
  double w;  /* state */
} zOscKura;

/* (static)
 * _zOscInitKura
 * - initialize Kuramoto's oscillator.
 */
void _zOscInitKura(void *prp, double p0, double w0)
{
  ((zOscKura *)prp)->p = p0;
  ((zOscKura *)prp)->w = w0;
}

/* (static)
 * _zOscUpdateKura
 * - update the internal state of Kuramoto's oscillator.
 */
double _zOscUpdateKura(void *prp, double u, double dt)
{
  zOscKura *kura;

  kura = prp;
  kura->w = kura->wn + kura->ke * sin( u - kura->p - kura->po );
  return kura->p += kura->w * dt;
}

/* (static)
 * _zOscOutputKura
 * - output a value of Kuramoto's oscillator.
 */
double _zOscOutputKura(void *prp)
{
  return ((zOscKura *)prp)->p;
}

static zOscCom _z_osc_com_kura = {
  _zOscInitKura,
  _zOscUpdateKura,
  _zOscOutputKura,
};

/* zOscCreateKura
 * - create Kuramoto's oscillator.
 */
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
