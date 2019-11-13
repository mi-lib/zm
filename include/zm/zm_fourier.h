/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_fourier - Fourier series.
 */

#ifndef __ZM_FOURIER_H__
#define __ZM_FOURIER_H__

#include <zm/zm_misc.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zFourier
 * (finite) Fourier series class
 * ********************************************************** */

typedef struct{
  int n;
  double *sc;
  double *cc;
} zFourier;

#define zFourierSetSinCoeffNC(f,i,c) ( (f)->sc[i] = (c) )
#define zFourierSetCosCoeffNC(f,i,c) ( (f)->cc[i] = (c) )
#define zFourierSetSinCoeff(f,i,c) ( (i) >= 0 && (i) < (f)->n ? zFourierSetSinCoeffNC( f, i, (i) == 0 ? 0 : (c) ) : 0)
#define zFourierSetCosCoeff(f,i,c) ( (i) >= 0 && (i) < (f)->n ? zFourierSetCosCoeffNC( f, i, c ) : 0 )

/*! \brief allocate internal memory for Fourier series. */
__EXPORT zFourier *zFourierAlloc(zFourier *f, int n);

/*! \brief free internal memory for Fourier series. */
__EXPORT void zFourierFree(zFourier *f);

/*! \brief value of a Fourier series. */
__EXPORT double zFourierVal(zFourier *f, double t);

/*! \brief velocity of a Fourier series. */
__EXPORT double zFourierVel(zFourier *f, double t);

/*! \brief acceleration of a Fourier series. */
__EXPORT double zFourierAcc(zFourier *f, double t);

__END_DECLS

#endif /* __ZM_FOURIER_H__ */
