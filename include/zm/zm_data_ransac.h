/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_data_ransac.h
 * \brief data analysis: random sampling consensus.
 * \author Zhidao
 */

#ifndef __ZM_DATA_RANSAC_H__
#define __ZM_DATA_RANSAC_H__

/* NOTE: never include this header file in user programs. */

#include <zm/zm_vec.h>

__BEGIN_DECLS

/*! \brief RANSAC: random sampling consensus.
 *
 * zRANSAC() identifies a model parameter \a q from the given samples
 * \a sample based on the model-fitting function \a fit_fp and the
 * error function \a error_fp. \a util is a utility pointer for
 * computations of the functions.
 * \a ns is the number of samples for a one-shot model guessing.
 * \a nt is the number of trials to find a candidate of the best parameter.
 * \a th is a threshold to distinguish inliers and outliers with respect
 * to the current guess.
 * \return
 * zRANSAC() returns the null pointer if it fails to allocate memory for
 * internal vector computations. Otherwise, a pointer \a q is returned.
 */
__EXPORT zVec zRANSAC(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,zVecList*,void*), double (* error_fp)(zVec,zVec,void*), void *util, int ns, int nt, double th);

__END_DECLS

#endif /* __ZM_DATA_RANSAC_H__ */
