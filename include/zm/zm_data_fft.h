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

#include <zm/zm_complex.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup Fast Fourier Transformation.
 * \{ *//* ************************************************** */

/*! \brief Fast Fourier Transformation.
 *
 * zFFT() is an implementation of the fast Fourier transformation
 * proposed by J.W.Cooley and J.W.Tukey(1965).
 * \a data is an array of sampled data with size \a size.
 * The result will be put into \a res as a series of complex
 * values. \a res[i] is the amplitude of the components with
 * frequency i/\a n[Hz]. The real part is the amplitude of
 * cosine components, while the imagenary part to that of sine
 * components.
 *
 * (Acknowledgment)
 * This implementation is based on the code by Dr. Takuya Ooura
 * (in http://www.kurims.kyoto-u.ac.jp/~ooura/fftman/index.html).
 *
 * \retval the false value if failing to internaly allocate a
 * working memory. Otherwise, the true value.
 */
__ZM_EXPORT bool zFFT(double data[], size_t n, zComplex res[]);

/*! \brief inverse fast Fourier transformation.
 *
 * zFFTInv() is an implementation of the inverse fast Fourier
 * transformation. \a data is an array of amplitudes with size \a n.
 * The result time series will be put into \a res.
 * Quantized time step between each component of \a res is 1/\a n.
 *
 * (Acknowledgement)
 * This implementation is based on the code by Dr. Takuya Ooura
 * (in http://www.kurims.kyoto-u.ac.jp/~ooura/fftman/index.html).
 *
 * \retval the false value if failing to internaly allocate
 * a working memory. Otherwise, the true value.
 */
__ZM_EXPORT bool zFFTInv(zComplex data[], size_t n, double res[]);

__END_DECLS

#endif /* __ZM_DATA_FFT_H__ */
