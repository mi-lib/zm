/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_data_fft.h
 * \brief data analysis: fast Fourier transformation.
 * \author Zhidao
 */

#ifndef __ZM_DATA_FFT_H__
#define __ZM_DATA_FFT_H__

/* NOTE: never include this header file in user programs. */

#include <zm/zm_cmat.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup Fast Fourier Transformation.
 * \{ *//* ************************************************** */

/*! \brief fast Fourier transformation.
 *
 * zFFT() is an implementation of the fast Fourier transformation proposed by J. W. Cooley and J. W. Tukey
 * (1965). \a data is a vector of sampled data.
 * The result spectrum will be put into a complex vector \a result. \a result[\a i] is the amplitude of the
 * components with frequency \a i/\a n[Hz]. The real part is the amplitude of cosine components, while the
 * imagenary part to that of sine components.
 *
 * zFFTInv() is an implementation of the inverse fast Fourier transformation. \a data is a complex spectrum
 * vector. The result time series will be put into a vector \a result.
 * Quantized time step between each component of \a result is 1/\a n.
 *
 * zFFT2() is an implementation of the two-dimensional fast Fourier transformation based on the Cooley-
 * Tukey algorithm. \a data is a matrix of sampled data. The result two-dimensional spectrum will be put
 * into a complex matrix \a result.
 *
 * zFFT2Inv() is an implementation of the two-dimensional inverse fast Fourier transformation. \a data is
 * a complex spectrum matrix. The result two-dimensional spacial series will be put into a matrix \a result.
 *
 * (Acknowledgment)
 * This implementation is based on the code by Dr. Takuya Ooura
 * (in http://www.kurims.kyoto-u.ac.jp/~ooura/fftman/index.html).
 * \return
 * These functions return the false value if the sizes of \a data and \a result are not equal or they fail
 * to allocate internal work space. Otherwise, the true value is returned.
 */
__ZM_EXPORT bool zFFT(zVec data, zCVec result);
__ZM_EXPORT bool zFFTInv(zCVec data, zVec result);
__ZM_EXPORT bool zFFT2(zMat data, zCMat result);
__ZM_EXPORT bool zFFT2Inv(zCMat data, zMat result);

__END_DECLS

#endif /* __ZM_DATA_FFT_H__ */
