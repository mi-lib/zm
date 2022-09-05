/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_data_ransac.h
 * \brief data analysis: random sample consensus.
 * \author Zhidao
 */

#ifndef __ZM_DATA_RANSAC_H__
#define __ZM_DATA_RANSAC_H__

/* NOTE: never include this header file in user programs. */

#include <zm/zm_vec.h>

__BEGIN_DECLS

/*! \brief RANSAC: random sample consensus.
 *
 * zRANSAC() identifies a model parameter \a q from the given samples
 * \a sample based on the model-fitting function \a fit_fp and the
 * error function \a error_fp. \a util is a utility pointer for
 * computations of the functions.
 * \a fit_fp( q, s, util ) should be a function that fits \a q to
 * the samples \a s given as a list of vectors.
 * \a error_fp( q, s, util ) should return the error between the true
 * value and an estimate value from the model \a q and a sample \a s.
 * \a ns is the number of samples for a one-shot model guessing.
 * \a nt is the number of trials to find a candidate of the best parameter.
 * \a th is a threshold to distinguish inliers and outliers with respect
 * to the current guess.
 *
 * zRANSACAuto() automatically sets parameters for RANSAC from the
 * rate of outliers and the level of noises.
 * \return
 * zRANSAC() and zRANSACAuto() return the null pointer if they fail to
 * allocate memory for internal vector computations. Otherwise, a pointer
 * \a q is returned.
 */
__EXPORT zVec zRANSAC(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,zVecList*,void*), double (* error_fp)(zVec,zVec,void*), void *util, uint ns, uint nt, double th);
__EXPORT zVec zRANSACAuto(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,zVecList*,void*), double (* error_fp)(zVec,zVec,void*), void *util, double r, double nl);

__END_DECLS

#endif /* __ZM_DATA_RANSAC_H__ */
