/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_fourier - Fourier series.
 */

#include <zm/zm_fourier.h>

/* allocate internal memory for Fourier series. */
zFourier *zFourierAlloc(zFourier *f, int n)
{
  f->n = n;
  f->sc = zAlloc( double, f->n );
  f->cc = zAlloc( double, f->n );
  if( !f->sc || !f->cc ){
    zFourierFree( f );
    return NULL;
  }
  return f;
}

/* free internal memory for Fourier series. */
void zFourierFree(zFourier *f)
{
  zFree( f->sc );
  zFree( f->cc );
  f->n = 0;
}

/* value of a Fourier series. */
double zFourierVal(zFourier *f, double t)
{
  int i;
  double s, c, val = 0;

  for( i=0; i<f->n; i++ ){
    zSinCos( i*t, &s, &c );
    val += f->sc[i]*s + f->cc[i]*c;
  }
  return val;
}

/* velocity of a Fourier series. */
double zFourierVel(zFourier *f, double t)
{
  int i;
  double s, c, val = 0;

  for( i=0; i<f->n; i++ ){
    zSinCos( i*t, &s, &c );
    val += i * ( f->sc[i]*c - f->cc[i]*s );
  }
  return val;
}

/* acceleration of a Fourier series. */
double zFourierAcc(zFourier *f, double t)
{
  int i;
  double s, c, val = 0;

  for( i=0; i<f->n; i++ ){
    zSinCos( i*t, &s, &c );
    val -= i*i * ( f->sc[i]*s + f->cc[i]*c );
  }
  return val;
}
