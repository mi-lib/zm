/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_mva_gmm.h
 * \brief multivariate analysis analysis: Gaussian mixture model
 * \author Zhidao
 */

#ifndef __ZM_MVA_GMM__
#define __ZM_MVA_GMM__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/*! \brief Gaussian model unit class
 *//* ******************************************************* */
typedef struct{
  zVec core;   /*!< core vector */
  zMat cov;    /*!< variance-covariance matrix */
  double weight;
  /*! \cond */
  zMat _cov_inv; /* inverse of covariance matrix */
  zVec _ci_e;    /* error multiplied by inverse of covariance */
  double _cov_det; /* determinant of covariance */
  /*! \endcond */
} zGMMUnit;

/*! \brief initialize a unit Gaussian model. */
__ZM_EXPORT zGMMUnit *zGMMUnitInit(zGMMUnit *gu);

/*! \brief allocate internal vectors and matrices of a unit Gaussian model. */
__ZM_EXPORT zGMMUnit *zGMMUnitAlloc(zGMMUnit *gu, int core_size, int error_size);

/*! \brief free internal vectors and matrices of a unit Gaussian model. */
__ZM_EXPORT void zGMMUnitFree(zGMMUnit *gu);

/* ********************************************************** */
/*! \brief Gaussian mixture model class
 *//* ******************************************************* */
zListClass( zGMMList, zGMMListCell, zGMMUnit );
typedef struct{
  zGMMList glist; /*!< list of Gaussian models */
  zVecClusterMethod method; /*!< methods for clustering */
  double log_likelihood; /*!< log-likelihood */
} zGMM;

/*! \brief initialize a Gaussian mixture model */
__ZM_EXPORT zGMM *zGMMInit(zGMM *gmm, int k, int size);

#define zGMMErrorF(gmm,v1,v2)          zVecClusterMethodErrorF( &(gmm)->method, v1, v2 )
#define zGMMDistF(gmm,v1,v2)           zVecClusterMethodDistF( &(gmm)->method, v1, v2 )
#define zGMMCoreF(gmm,vl,c)            zVecClusterMethodCoreF( &(gmm)->method, vl, c )
#define zGMMLoadedMeanF(gmm,vl,l,nk,m) zVecClusterMethodLoadedMeanF( &(gmm)->method, vl, l, nk, m )

__ZM_EXPORT bool zGMMSetErrorFunc(zGMM *gmm, int error_size, zVec (* error_fp)(const zVecClusterMethod*,const zVec,const zVec,void*,zVec), void *util);
__ZM_EXPORT bool zGMMSetDistFunc(zGMM *gmm, double (* dist_fp)(const zVecClusterMethod*,const zVec,const zVec,void*), void *util);
__ZM_EXPORT bool zGMMSetCoreFunc(zGMM *gmm, int core_size, zVec (* core_fp)(const zVecClusterMethod*,const zVecAddrList*,void*,zVec), void *util);
__ZM_EXPORT bool zGMMSetLoadedMeanFunc(zGMM *gmm, zVec (* lm_fp)(const zVecClusterMethod*,const zVecAddrList*,const double[],double,void*,zVec), void *util);

/*! \brief set functions to cluster data based on the linear-sum for Gaussian mixture model. */
__ZM_EXPORT bool zGMMSetMethodLS(zGMM *gmm, void *util);

/*! \brief destroy a Gaussian mixture model */
__ZM_EXPORT void zGMMDestroy(zGMM *gmm);

/*! \brief create a Gaussian mixture model based on EM algorithm */
__ZM_EXPORT zGMM *zGMMCreateEM(zGMM *gmm, const zVecList *points);

/*! \brief Akaike's Information Criterion. */
__ZM_EXPORT double zGMMAIC(const zGMM *gmm);

/*! \brief Bayesian Information Criterion. */
__ZM_EXPORT double zGMMBIC(const zGMM *gmm, const zVecAddrList *sample);

__END_DECLS

#endif /* __ZM_MVA_GMM__ */
