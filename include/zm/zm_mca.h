/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_mca.h
 * \brief multiple classification analysis.
 * \author Zhidao
 */

#ifndef __ZM_MCA_H__
#define __ZM_MCA_H__

#include <zm/zm_eig.h>

__BEGIN_DECLS

/*! \brief sum up all vectors in a list.
 */
__EXPORT zVec zVecListSum(zVecList *list, zVec sum);

/*! \brief average of all vectors in a list.
 */
__EXPORT zVec zVecListAve(zVecList *list, zVec ave);

/*! \brief variance of all vectors in a list.
 */
__EXPORT double zVecListVar(zVecList *list, zVec ave);

/*! \brief average and variance of all vectors in a list.
 */
__EXPORT double zVecListAveVar(zVecList *list, zVec ave);

/*! \brief variance-covariance matrix of all vectors in a list.
 */
__EXPORT zMat zVecListCov(zVecList *list, zVec ave, zMat cov);

/*! \brief average and variance-covariance matrix of all vectors in a list.
 */
__EXPORT zMat zVecListAveCov(zVecList *list, zVec ave, zMat cov);

/*! \brief principal component analysis.
 */
__EXPORT int zPCA(zVecList *points, double cr, zVec ave, zVec score, zMat loading);

__END_DECLS

#include <zm/zm_mca_cluster.h>
#include <zm/zm_mca_gmm.h>

#endif /* __ZM_MCA_H__ */
