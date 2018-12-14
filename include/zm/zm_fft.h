/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_fft.h
 * \brief Fast Fourier Transform.
 * \author Zhidao
 */

#ifndef __Z_FFT_H__
#define __Z_FFT_H__

#include <zm/zm_complex.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup Fast Fourier Transform.
 * \{ *//* ************************************************** */

/*! \brief Fast Fourier Transform.
 *
 * zFFT() is an implementation of the Fast Fourier Transform
 * proposed by J.W.Cooley and J.W.Tukey(1965).
 * \a data is an array of sampled data with size \a n.
 * The result will be put into \a res as a series of complex
 * values. \a res[i] is the amplitude of the components with
 * frequency i/\a n[Hz]. The real part is the amplitude of
 * cosine components, while the imagenary part to that of
 * sine components.
 *
 * \retval the false value if failing to internaly allocate
 * a working memory. Otherwise, the true value.
 * [ACKNOWLEDGEMENT]
 * This implementation is based on the code by Mr.Takuya Ooura
 * (in http://www.kurims.kyoto-u.ac.jp/~ooura/fftman/index.html).
 */
__EXPORT bool zFFT(double data[], int n, zComplex res[]);

/*! \brief inverse Fast Fourier Transform.
 *
 * zFFTInv() is an implementation of the inverse Fast Fourier
 * Transform. \a data is an array of amplitudes with size \a n.
 * The result time series will be put into \a res.
 * Quantized time step between each component of \a res is 1/\a n.
 *
 * \retval the false value if failing to internaly allocate
 * a working memory. Otherwise, the true value.
 * [ACKNOWLEDGEMENT]
 * This implementation is based on the code by Mr.Takuya Ooura
 * (in http://www.kurims.kyoto-u.ac.jp/~ooura/fftman/index.html).
 */
__EXPORT bool zFFTInv(zComplex data[], int n, zComplex res[]);

__END_DECLS

#endif /* __Z_FFT_H__ */
