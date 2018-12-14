/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_mca_cluster.h
 * \brief multiple classification analysis : clustering.
 * \author Zhidao
 */

#ifndef __ZM_MCA_CLUSTER_H__
#define __ZM_MCA_CLUSTER_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/*! \brief vecter cluster class.
 *//* ******************************************************* */
typedef struct{
  zVecList vl; /*!< \brief list of vectors */
  zVec mean;   /*!< \brief mean vector of the cluster */
  double var;  /*!< \brief variance */
} zCluster;

/*! \brief create a vector cluster */
__EXPORT zCluster *zClusterCreate(zCluster *c, int meansize);

/*! \brief destroy a vector cluster */
__EXPORT void zClusterDestroy(zCluster *c);

/*! \brief output a vector cluster to a file */
__EXPORT void zClusterFWrite(FILE *fp, zCluster *c);
__EXPORT void zClusterDataFWrite(FILE *fp, zCluster *c);

/* ********************************************************** */
/*! \brief methods for mean and error computation
 *//* ******************************************************* */

typedef struct{
  /*! \cond */
  int _meansize; /* size of the mean vector */
  int _errorsize; /* size of the error vector */
  zVec (* _mean_fp)(zVecList*,void*,zVec); /* mean function */
  zVec (* _error_fp)(zVec,zVec,void*,zVec); /* error function */
  zVec _err; /* error vector */

  zVec (* _mean_l_fp)(zVecList*,double[],double,void*,zVec); /* loaded mean function */
  zVec _err_l; /* loaded error vector */
  /*! \endcond */
} zClusterMethod;

__EXPORT void zClusterMethodInit(zClusterMethod *met);

__EXPORT zClusterMethod *zClusterMethodCreate(zClusterMethod *met, int meansize, zVec (* mean_fp)(zVecList*,void*,zVec), int errorsize, zVec (* error_fp)(zVec,zVec,void*,zVec));
__EXPORT zClusterMethod *zClusterMethodLoadedCreate(zClusterMethod *met, zVec (* mean_l_fp)(zVecList*,double[],double,void*,zVec));
__EXPORT void zClusterMethodDestroy(zClusterMethod *met);

/* ********************************************************** */
/*! \brief multiple vecter clusters class.
 *//* ******************************************************* */

zListClass( zClusterList, zClusterListCell, zCluster );

typedef struct{
  zClusterList cl; /*!< a list of clusters */
  zClusterMethod met; /*!< methods for mean and error computation */
} zMCluster;

/*! \brief initialize multiple vector clusters */
__EXPORT zMCluster *zMClusterInit(zMCluster *mc, int meansize, zVec (* mean_fp)(zVecList*,void*,zVec), int errorsize, zVec (* error_fp)(zVec,zVec,void*,zVec));

/*! \brief allocate a set of vector clusters */
__EXPORT zMCluster *zMClusterAlloc(zMCluster *mc, int n);

/*! \brief destroy a set of vector clusters */
__EXPORT void zMClusterDestroy(zMCluster *mc);

/*! \brief output vectors in a set of vector clusters */
__EXPORT void zMClusterFWrite(FILE *fp, zMCluster *mc);
__EXPORT void zMClusterDataFWrite(FILE *fp[], zMCluster *mc);
__EXPORT void zMClusterMeanFWrite(FILE *fp[], zMCluster *mc);

/* ********************************************************** */
/*! \brief clustering based on K-means family
 *//* ******************************************************* */

/*! \brief clustering of vectors by K-means */
__EXPORT int zMClusterKMeans(zMCluster *mc, zVecList *points, int k, void *mean_util, void *err_util);

/*! \brief clustering of vectors by X-means */
__EXPORT int zMClusterXMeans(zMCluster *mc, zVecList *points, void *mean_util, void *err_util);

/*! \brief clustering of vectors by X-means based on BIC */
__EXPORT int zMClusterXMeansBIC(zMCluster *mc, zVecList *points, void *mean_util, void *err_util);

__END_DECLS

#endif /* __ZM_MCA_CLUSTER_H__ */
