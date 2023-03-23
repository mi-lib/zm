/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_mca_gmm.h
 * \brief multiple classification analysis: Gaussian mixture model
 * \author Zhidao
 */

#ifndef __ZM_MCA_GMM__
#define __ZM_MCA_GMM__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/*! \brief unit Gaussian function class
 *//* ******************************************************* */
typedef struct{
  bool active; /*!< flag to activate */
  zVec mean; /*!< mean vector */
  zMat cov;  /*!< variance-covariance matrix */
  double weight;
  /*! \cond */
  zMat _cov_inv; /* inverse of covariance matrix */
  zVec _ci_e;    /* error multiplied by inverse of covariance */
  double _cov_det; /* determinant of covariance */
  /*! \endcond */
} zGMMUnit;

/*! \brief allocate internal vectors and matrices of a unit Gaussian function */
__EXPORT zGMMUnit *zGMMUnitAlloc(zGMMUnit *gu, int meansize, int errorsize);

/*! \brief free internal vectors and matrices of a unit Gaussian function */
__EXPORT void zGMMUnitFree(zGMMUnit *gu);

/* ********************************************************** */
/*! \brief Gaussian mixture model class
 *//* ******************************************************* */
zListClass( zGMMList, zGMMListCell, zGMMUnit );
typedef struct{
  zGMMList gl; /*!< list of Gaussian functions */
  zClusterMethod method; /*!< methods for clustering */
  double log_likelihood; /*!< log-likelihood */
} zGMM;

/*! \brief initialize a Gaussian mixture model */
__EXPORT zGMM *zGMMInit(zGMM *gmm, int k, int meansize, zVec (* mean_fp)(zClusterMethod*,zVecList*,void*,zVec), void *mean_util, zVec (* l_mean_fp)(zClusterMethod*,zVecList*,double[],double,void*,zVec), void *l_mean_util, int errorsize, zVec (* error_fp)(zClusterMethod*,zVec,zVec,void*,zVec), void *error_util, double (* dist_fp)(zClusterMethod*,zVec,zVec,void*), void *dist_util);

/*! \brief destroy a Gaussian mixture model */
__EXPORT void zGMMDestroy(zGMM *gmm);

/*! \brief create a Gaussian mixture model based on EM algorithm */
__EXPORT zGMM *zGMMCreateEM(zGMM *gmm, zVecList *points, int k);

__END_DECLS

#endif /* __ZM_MCA_GMM__ */
