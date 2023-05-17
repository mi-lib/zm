/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_data.h
 * \brief data analysis.
 * \author Zhidao
 */

#ifndef __ZM_DATA_H__
#define __ZM_DATA_H__

#include <zm/zm_misc.h>

__BEGIN_DECLS

/*! \brief pick up peaks a data sequence. */
__ZM_EXPORT zIndex zDataPeak(double src[], size_t n, uint w);

/*! \brief smooth a data sequence based on Savitzky-Golay's method */
__ZM_EXPORT bool zDataSmoothSG(double src[], size_t n, size_t w, uint dim, double dest[]);

/*! \brief pick up peaks of a smoothed data sequence based on Savitzky-Golay's method. */
__ZM_EXPORT zIndex zDataPeakSG(double src[], size_t n, size_t w, uint dim);

/*! \brief sort an array of double-precision floating-point vector in the descent order. */
__ZM_EXPORT double *zDataSort(double data[], size_t n);

/*! \brief sort an integer vector in the descent order of the corresponding samples. */
__ZM_EXPORT zIndex zDataSortIndex(double data[], size_t n, zIndex index);

/*! \brief sort an integer vector in the descent order of absolute values of the corresponding samples. */
__ZM_EXPORT zIndex zDataSortAbsIndex(double data[], size_t n, zIndex index);

__END_DECLS

#include <zm/zm_data_fft.h>
#include <zm/zm_data_ransac.h>

#endif /* __ZM_DATA_H__ */
