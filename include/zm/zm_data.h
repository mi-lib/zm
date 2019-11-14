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

/*! \brief smooth a data sequence based on Savitzky-Golay's method */
__EXPORT bool zDataSmoothSG(double src[], size_t n, size_t w, int dim, double dest[]);

/*! \brief smooth a data sequence based on Savitzky-Golay's method and pick up peaks. */
__EXPORT zIndex zDataPeakSG(double src[], size_t n, int w, int dim);

__END_DECLS

#include <zm/zm_data_fft.h>

#endif /* __ZM_DATA_H__ */
