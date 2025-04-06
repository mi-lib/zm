/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_mva_ransac.h
 * \brief  multivariate analysis analysis: random sample consensus.
 * \author Zhidao
 */

#ifndef __ZM_MVA_RANSAC_H__
#define __ZM_MVA_RANSAC_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief RANSAC: random sample consensus.
 *
 * zRANSAC() identifies a model parameter \a q from the given samples \a sample based on the model-fitting
 * function \a fit_fp and the error function \a error_fp. \a util is a utility pointer for computations of
 * the functions.
 * \a fit_fp( q, s, util ) should be a function that fits \a q to the samples \a s given as a list of vectors.
 * \a error_fp( q, s, util ) should return the error between the true value and an estimate value from the
 * model \a q and a sample \a s.
 * \a ns is the number of samples for a one-shot model guessing.
 * \a nt is the number of trials to find a candidate of the best parameter.
 * \a th is a threshold to distinguish inliers and outliers with respect to the current guess.
 *
 * zRANSACSaveInlier() works in the same way with zRANSAC() except that it saves the candidate samples of
 * inliers to a list of vectors \a inlier from \a sample.
 *
 * zRANSACAuto() and zRANSACSaveInlierAuto() automatically set parameters for RANSAC from the rate of
 * outliers \a rate and the level of noises \a nl. Their arguments have the same meanings with zRANSAC()
 * and zRANSACSaveInlier(), respectively.
 * \return
 * zRANSAC() and zRANSACAuto() return the null pointer if they fail to allocate memory for internal vector
 * computations. Otherwise, a pointer \a q is returned.
 */
__ZM_EXPORT zVec zRANSAC(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,const zVecList*,void*), double (* error_fp)(const zVec,const zVec,void*), void *util, int ns, int nt, double th);
__ZM_EXPORT zVec zRANSACSaveInlier(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,const zVecList*,void*), double (* error_fp)(const zVec,const zVec,void*), void *util, int ns, int nt, double th, zVecList *inlier);
__ZM_EXPORT zVec zRANSACAuto(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,const zVecList*,void*), double (* error_fp)(const zVec,const zVec,void*), void *util, double rate, double nl);
__ZM_EXPORT zVec zRANSACSaveInlierAuto(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,const zVecList*,void*), double (* error_fp)(const zVec,const zVec,void*), void *util, double rate, double nl, zVecList *inlier);

__END_DECLS

#endif /* __ZM_MVA_RANSAC_H__ */
