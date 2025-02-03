/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_fourier - Fourier series.
 */

#ifndef __ZM_FOURIER_H__
#define __ZM_FOURIER_H__

#include <zm/zm_misc.h>

__BEGIN_DECLS

/*! \struct zFourier
 * \brief finite Fourier series class
 */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zFourier ){
  int n;      /*!< \brief number of terms */
  double *sc; /*!< \brief coefficients on sine terms */
  double *cc; /*!< \brief coefficients on cosine terms */
};

/*! \brief set a coefficient on a sine term without checking the index. */
#define zFourierSetSinCoeffNC(f,i,c) ( (f)->sc[i] = (c) )
/*! \brief set a coefficient on a cosine term without checking the index. */
#define zFourierSetCosCoeffNC(f,i,c) ( (f)->cc[i] = (c) )

/*! \brief set a coefficient on a sine term. */
#define zFourierSetSinCoeff(f,i,c) ( (i) >= 0 && (i) < (f)->n ? zFourierSetSinCoeffNC( f, i, (i) == 0 ? 0 : (c) ) : 0)
/*! \brief set a coefficient on a cosine term. */
#define zFourierSetCosCoeff(f,i,c) ( (i) >= 0 && (i) < (f)->n ? zFourierSetCosCoeffNC( f, i, c ) : 0 )

/*! \brief allocate internal memory for Fourier series. */
__ZM_EXPORT zFourier *zFourierAlloc(zFourier *f, int n);

/*! \brief free internal memory for Fourier series. */
__ZM_EXPORT void zFourierFree(zFourier *f);

/*! \brief value of a Fourier series. */
__ZM_EXPORT double zFourierVal(zFourier *f, double t);

/*! \brief velocity of a Fourier series. */
__ZM_EXPORT double zFourierVel(zFourier *f, double t);

/*! \brief acceleration of a Fourier series. */
__ZM_EXPORT double zFourierAcc(zFourier *f, double t);

/*! \brief print out coefficients of a Fourier series. */
__ZM_EXPORT void zFourierFPrint(FILE *fp, zFourier *fourier);

__END_DECLS

#endif /* __ZM_FOURIER_H__ */
