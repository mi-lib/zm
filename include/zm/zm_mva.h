/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_mva.h
 * \brief multivariate analysis analysis.
 * \author Zhidao
 */

#ifndef __ZM_MVA_H__
#define __ZM_MVA_H__

#include <zm/zm_data.h>
#include <zm/zm_mat_eig.h>

__BEGIN_DECLS

/*! \brief minimum and maximum of all vectors in a list.
 *
 * zVecListMinMax() finds the minimum and maximum values of each component of vectors in a list \a list
 * and put them into \a min and \a max, respectively.
 * In other words, vectors \a min and \a max are composed from the minimum and maximum values of each
 * component.
 * \return
 * zVecListMinMax() returns the true value if the sizes of \a min, \a max and vectors in \a list have
 * the same size. Otherwise, it returns the false value.
 */
__ZM_EXPORT bool zVecListMinMax(const zVecList *list, zVec min, zVec max);

/*! \brief sum up all vectors in a list.
 *
 * zVecListSum() sums up all vectors of a list \a list to \a sum.
 * \return
 * zVecListSum() returns a pointer \a sum.
 */
__ZM_EXPORT zVec zVecListSum(const zVecList *list, zVec sum);

/*! \brief mean of all vectors in a list.
 *
 * zVecListMean() computes the mean vector \a mean of a list \a list.
 * \return
 * zVecListMean() returns a pointer \a mean.
 */
__ZM_EXPORT zVec zVecListMean(const zVecList *list, zVec mean);

/*! \brief variance of all vectors in a list.
 *
 * zVecListVar() computes the variance of all vectors of a list \a list with respect to a given mean
 * vector \a mean.
 * \return
 * zVecListVar() returns the computed variance.
 */
__ZM_EXPORT double zVecListVar(const zVecList *list, zVec mean);

/*! \brief mean and variance of all vectors in a list.
 *
 * zVecListMeanVar() simultaneously computes the mean vector and the variance of a list \a list. The
 * former is stored in \a mean and the latter is returned by the function.
 * \return
 * zVecListMeanVar() returns the computed variance.
 */
__ZM_EXPORT double zVecListMeanVar(const zVecList *list, zVec mean);

/*! \brief variance-covariance matrix of all vectors in a list.
 *
 * zVecListCov() computes the variance-covariance matrix of all vectors of a list \a list with respect
 * to a given mean vector \a mean. The result is stored in \a cov.
 * \return
 * zVecListCov() returns a pointer \a cov.
 */
__ZM_EXPORT zMat zVecListCov(const zVecList *list, zVec mean, zMat cov);

/*! \brief mean and variance-covariance matrix of all vectors in a list.
 *
 * zVecListMeanCov() simultaneously computes the mean vector and the variance-covariance matrix of
 * a list \a list. The former is stored in \a mean and the latter in \a cov.
 * \return
 * zVecListMeanCov() returns a pointer \a cov.
 */
__ZM_EXPORT zMat zVecListMeanCov(const zVecList *list, zVec mean, zMat cov);

/*! \brief principal component analysis.
 *
 * zVecListPCA() computes the mean vector, the score and the loading of a give set of vector points
 * of a list \a points through the principal component analysis, and store them into \a mean, \a score,
 * and \a loading, respectively.
 * \a cr is the threshold contribution ratio. If the score of a component is less than it, the component
 * is ignored.
 * \return
 * zVecListPCA() returns the number of significant (i.e., unignored) components.
 */
__ZM_EXPORT int zVecListPCA(const zVecList *points, double cr, zVec mean, zVec score, zMat loading);

/*! \brief generate vectors from normal distribution.
 *
 * zVecListGenRandND() randomly generates vectors that follow normal distribution defined by a mean
 * vector \a mean and a variance-covariance matrix \a cov. \a n is the number of vectors to be generated.
 * The results are stored in a list \a vl.
 * This function fails if the sizes of \a mean and \a cov are inconsistent.
 * \return
 * zVecListGenRandND() returns the number of actually generated vectors, namely, the size of \a vl.
 * If it internally fails to allocate memory, a smaller number than \a n is returned.
 */
__ZM_EXPORT int zVecListGenRandND(zVecList *vl, int n, zVec mean, zMat cov);

__END_DECLS

#include <zm/zm_mva_ransac.h>
#include <zm/zm_mva_cluster.h>
#include <zm/zm_mva_gmm.h>

#endif /* __ZM_MVA_H__ */
