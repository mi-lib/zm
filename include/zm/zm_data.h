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
__ZM_EXPORT zIndex zDataPeak(double src[], int size, int w);

/*! \brief sort an array of double-precision floating-point vector in the descent order. */
__ZM_EXPORT double *zDataSort(double data[], int size);

/*! \brief sort an integer vector in the descent order of the corresponding samples. */
__ZM_EXPORT zIndex zDataSortIndex(double data[], int size, zIndex index);

/*! \brief sort an integer vector in the descent order of absolute values of the corresponding samples. */
__ZM_EXPORT zIndex zDataSortAbsIndex(double data[], int size, zIndex index);

/*! \brief select data of an array of double-precision floating-point values at the designated position in the ascent order. */
__ZM_EXPORT double zDataSelect(double data[], int size, int select_order);

/* median of an array of double-precision floating-point values. */
__ZM_EXPORT double zDataMedian(double data[], int size);

/*! \brief smooth a data sequence based on Savitzky-Golay's method */
__ZM_EXPORT bool zDataSmoothSG(double src[], int size, int w, int dim, double dest[]);

/*! \brief pick up peaks of a smoothed data sequence based on Savitzky-Golay's method. */
__ZM_EXPORT zIndex zDataPeakSG(double src[], int size, int w, int dim);

__END_DECLS

#include <zm/zm_data_fft.h>

#endif /* __ZM_DATA_H__ */
