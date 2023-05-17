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
  zVecAddrList vlist; /*!< \brief list of vectors */
  zVec core;          /*!< \brief core vector of the cluster */
  double var;         /*!< \brief variance */
  /*! \cond */
  double *_sil;       /* silhouettes */
  /*! \endcond */
} zCluster;

#define zClusterSampleList(c) ( &(c)->vlist )
#define zClusterCore(c)       (c)->core
#define zClusterVar(c)        (c)->var

/*! \brief create a vector cluster */
__ZM_EXPORT zCluster *zClusterCreate(zCluster *c, int coresize);

/*! \brief destroy a vector cluster */
__ZM_EXPORT void zClusterDestroy(zCluster *c);

/*! \brief the maximum silhouette in a cluster. */
__ZM_EXPORT double zClusterMaxSilhouette(zCluster *c);

/*! \brief print a vector cluster to a file */
__ZM_EXPORT void zClusterFPrint(FILE *fp, zCluster *c);
__ZM_EXPORT void zClusterDataFPrint(FILE *fp, zCluster *c);

/* ********************************************************** */
/*! \brief methods for core and error computation
 *//* ******************************************************* */

ZDEF_STRUCT( zClusterMethod ){
  /*! \cond */
  zVec error; /* error vector */
  zVec (* error_fp)(zClusterMethod*,zVec,zVec,void*,zVec); /* error function */
  void *error_util;
  double (* dist_fp)(zClusterMethod*,zVec,zVec,void*); /* distance function */
  void *dist_util;
  int core_size; /* size of the core vector */
  zVec (* core_fp)(zClusterMethod*,zVecAddrList*,void*,zVec); /* a function to find core */
  void *core_util;
  /* for GMM */
  zVec (* lm_fp)(zClusterMethod*,zVecAddrList*,double[],double,void*,zVec); /* loaded mean function */
  void *lm_util;
  zVec _lerr; /* workspace for loaded error vector */
  /*! \endcond */
};

#define zClusterMethodErrorF(cm,v1,v2) (cm)->error_fp( cm, (v1), (v2), (cm)->error_util, (cm)->error )
#define zClusterMethodDistF(cm,v1,v2)  (cm)->dist_fp( cm, (v1), (v2), (cm)->dist_util )
#define zClusterMethodCoreF(cm,vl,c)   (cm)->core_fp( cm, (vl), (cm)->core_util, (c) )
#define zClusterMethodLoadedMeanF(cm,vl,l,nk,m) (cm)->lm_fp( cm, (vl), (l), (nk), (cm)->lm_util, (m) )

/*! \brief initialize methods for clustering. */
__ZM_EXPORT void zClusterMethodInit(zClusterMethod *method);

/*! \brief set a function to find core of a cluster. */
__ZM_EXPORT zClusterMethod *zClusterMethodSetCoreFunc(zClusterMethod *method, int size, zVec (* fp)(zClusterMethod*,zVecAddrList*,void*,zVec), void *util);
/*! \brief set a function to compute error of a sample in a cluster. */
__ZM_EXPORT zClusterMethod *zClusterMethodSetErrorFunc(zClusterMethod *method, int size, zVec (* fp)(zClusterMethod*,zVec,zVec,void*,zVec), void *util);
/*! \brief set a function to compute the distance between two samples in a cluster. */
__ZM_EXPORT zClusterMethod *zClusterMethodSetDistFunc(zClusterMethod *method, double (* fp)(zClusterMethod*,zVec,zVec,void*), void *util);
/*! \brief set a function to compute the loaded mean of samples (for GMM). */
__ZM_EXPORT zClusterMethod *zClusterMethodSetLoadedMeanFunc(zClusterMethod *method, zVec (* mean_l_fp)(zClusterMethod*,zVecAddrList*,double[],double,void*,zVec), void *util);

__ZM_EXPORT zClusterMethod *zClusterMethodCreate(zClusterMethod *method, int size);

/*! \brief copy methods for clustering. */
__ZM_EXPORT zClusterMethod *zClusterMethodCopy(zClusterMethod *src, zClusterMethod *dest);

/*! \brief destroy methods for clustering. */
__ZM_EXPORT void zClusterMethodDestroy(zClusterMethod *method);

/* ********************************************************** */
/*! \brief multiple vecter clusters class.
 *//* ******************************************************* */

zListClass( zClusterList, zClusterListCell, zCluster );

ZDEF_STRUCT( zMCluster ){
  zClusterList clist;    /*!< a list of clusters */
  zClusterMethod method; /*!< methods for core and error computation */
};

#define zMClusterClusterList(mc) ( &(mc)->clist )

/*! \brief initialize multiple vector clusters.
 *
 * zMClusterInit() intializes a multiple clusters \a mc.
 * \a coresize is the size of core vector, i.e., adjoint vector to samples, of a cluster.
 * \a errorsize is the size of error vector.
 * In typical cases where the core means the mean of samples, both \a coresize and
 * \a errorsize are equal to the dimension of samples.
 * \return
 * zMClusterInit() returns a pointer \a mc if it succeeds to initialize \a mc.
 * Otherwise, the null pointer is returned.
 */
__ZM_EXPORT zMCluster *zMClusterInit(zMCluster *mc, int size);

#define zMClusterSetErrorFunc(mc,se,fp,util) zClusterMethodSetErrorFunc( &(mc)->method, (se), (fp), (util) )
#define zMClusterSetDistFunc(mc,fp,util)     zClusterMethodSetDistFunc( &(mc)->method, (fp), (util) )
#define zMClusterSetCoreFunc(mc,sc,fp,util)  zClusterMethodSetCoreFunc( &(mc)->method, (sc), (fp), (util) )

#define zMClusterErrorF(mc,v1,v2) zClusterMethodErrorF( &(mc)->method, v1, v2 )
#define zMClusterDistF(mc,v1,v2)  zClusterMethodDistF( &(mc)->method, v1, v2 )
#define zMClusterCoreF(mc,vl,c)   zClusterMethodCoreF( &(mc)->method, vl, c )

/*! \brief allocate a set of vector clusters.
 *
 * zMClusterAlloc() allocates memory for multiple clusters \a mc.
 * \a n is the number of clusters.
 * \return
 * zMClusterAlloc() returns the pointer \a mc if it succeeds to allocate
 * memory. Otherwise, the null pointer is returned.
 */
__ZM_EXPORT zMCluster *zMClusterAlloc(zMCluster *mc, int n);

/*! \brief copy methods of a multiple cluster to another. */
#define zMClusterMethodCopy(src,dest) zClusterMethodCopy( &(src)->method, &(dest)->method )

/*! \brief move a multiple cluster to another. */
__ZM_EXPORT bool zMClusterMove(zMCluster *src, zMCluster *dest);

/*! \brief destroy a set of vector clusters. */
__ZM_EXPORT void zMClusterDestroy(zMCluster *mc);

/*! \brief evenness of clusters.
 *
 * zMClusterEvenness() returns the measure of evenness, namely, the ratio of
 * the maximum number of clustered samples over the minimum, of a multiple
 * cluster \a mc. \a mc is more even if it has as close value to 1.
 * \return
 * zMClusterEvenness() returns the measure of evenness.
 */
__ZM_EXPORT double zMClusterEvenness(zMCluster *mc);

/*! \brief print vectors in a set of vector clusters */
__ZM_EXPORT void zMClusterFPrint(FILE *fp, zMCluster *mc);
__ZM_EXPORT void zMClusterDataFPrint(FILE *fp[], zMCluster *mc);
__ZM_EXPORT void zMClusterCoreFPrint(FILE *fp[], zMCluster *mc);

/*! \brief print vectors in a set of clusters to files with a common basename. */
__ZM_EXPORT bool zMClusterDataPrintFile(zMCluster *mc, const char *basename);

/* ********************************************************** */
/*! \brief clustering based on K-means
 *//* ******************************************************* */

/*! \brief clustering of vectors by K-means.
 *
 * zMClusterKMeansKKZ() and zMClusterKMeans() make multiple clusters from a set of
 * samples based on K-means method. zMClusterKMeansKKZ() uses KKZ method for
 * initialization, while zMClusterKMeans() uses K-means++.
 * For both functions, the resulted clusters are stored in \a mc. \a points is the
 * list of pointers to the original samples to be clusters. \a k is the number of
 * clusters.
 * \return
 * zMClusterKMeansKKZ() and zMClusterKMeans() return the number of iterations that
 * is taken in the K-means method.
 */
__ZM_EXPORT int zMClusterKMeansKKZ(zMCluster *mc, zVecAddrList *points, int k);
__ZM_EXPORT int zMClusterKMeans(zMCluster *mc, zVecAddrList *points, int k);

/*! \brief clustering of vectors by K-medoids.
 *
 * zMClusterKMedoids() makes multiple clusters from a set of samples based on K-medoids
 * method. It initializes clusters by the same method with K-means++.
 * The resulted clusters are stored in \a mc. \a points is the list of pointers to
 * the original samples to be clusters. \a k is the number of clusters.
 * \return
 * zMClusterKMedoids() returns the number of iterations that is taken in the method.
 * \sa
 * zMClusterKMeansKKZ, zMClusterKMeans
 */
__ZM_EXPORT int zMClusterKMedoids(zMCluster *mc, zVecAddrList *points, int k);

/*! \brief compute the mean silhouette of a set of clusters. */
__ZM_EXPORT double zMClusterMeanSilhouette(zMCluster *mc);

/*! \brief print silhouettes of a set of vector clusters to files. */
__ZM_EXPORT bool zMClusterSilhouettePrintFile(zMCluster *mc, const char *basename);

/* ********************************************************** */
/*! \brief clustering based on X-means
 *//* ******************************************************* */

/*! \brief clustering of vectors by X-means based on hyperdensity of clusters. */
__ZM_EXPORT int zMClusterXMeansDensity(zMCluster *mc, zVecAddrList *points);

/*! \brief clustering of vectors by X-means based on BIC */
__ZM_EXPORT int zMClusterXMeansBIC(zMCluster *mc, zVecAddrList *points);

__END_DECLS

#endif /* __ZM_MCA_CLUSTER_H__ */
