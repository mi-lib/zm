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
__ZM_EXPORT zGMMUnit *zGMMUnitAlloc(zGMMUnit *gu, int coresize, int errorsize);

/*! \brief free internal vectors and matrices of a unit Gaussian model. */
__ZM_EXPORT void zGMMUnitFree(zGMMUnit *gu);

/* ********************************************************** */
/*! \brief Gaussian mixture model class
 *//* ******************************************************* */
zListClass( zGMMList, zGMMListCell, zGMMUnit );
typedef struct{
  zGMMList glist; /*!< list of Gaussian models */
  zClusterMethod method; /*!< methods for clustering */
  double log_likelihood; /*!< log-likelihood */
} zGMM;

/*! \brief initialize a Gaussian mixture model */
__ZM_EXPORT zGMM *zGMMInit(zGMM *gmm, int k, int size);

__ZM_EXPORT bool zGMMSetErrorFunc(zGMM *gmm, int size, zVec (* error_fp)(zClusterMethod*,zVec,zVec,void*,zVec), void *util);
__ZM_EXPORT bool zGMMSetDistFunc(zGMM *gmm, double (* dist_fp)(zClusterMethod*,zVec,zVec,void*), void *util);
__ZM_EXPORT bool zGMMSetCoreFunc(zGMM *gmm, int size, zVec (* core_fp)(zClusterMethod*,zVecAddrList*,void*,zVec), void *util);
__ZM_EXPORT bool zGMMSetLoadedMeanFunc(zGMM *gmm, zVec (* lm_fp)(zClusterMethod*,zVecAddrList*,double[],double,void*,zVec), void *util);

#define zGMMErrorF(gmm,v1,v2)          zClusterMethodErrorF( &(gmm)->method, v1, v2 )
#define zGMMDistF(gmm,v1,v2)           zClusterMethodDistF( &(gmm)->method, v1, v2 )
#define zGMMCoreF(gmm,vl,c)            zClusterMethodCoreF( &(gmm)->method, vl, c )
#define zGMMLoadedMeanF(gmm,vl,l,nk,m) zClusterMethodLoadedMeanF( &(gmm)->method, vl, l, nk, m )

/*! \brief destroy a Gaussian mixture model */
__ZM_EXPORT void zGMMDestroy(zGMM *gmm);

/*! \brief create a Gaussian mixture model based on EM algorithm */
__ZM_EXPORT zGMM *zGMMCreateEM(zGMM *gmm, zVecList *points);

/*! \brief Akaike's Information Criterion. */
__ZM_EXPORT double zGMMAIC(zGMM *gmm);

/*! \brief Bayesian Information Criterion. */
__ZM_EXPORT double zGMMBIC(zGMM *gmm, zVecAddrList *sample);

__END_DECLS

#endif /* __ZM_MCA_GMM__ */
