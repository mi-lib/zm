/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_mva_cluster.h
 * \brief multivariate analysis analysis : clustering.
 * \author Zhidao
 */

#ifndef __ZM_MVA_CLUSTER_H__
#define __ZM_MVA_CLUSTER_H__

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
} zVecCluster;

#define zVecClusterSampleList(c) ( &(c)->vlist )
#define zVecClusterCore(c)       (c)->core
#define zVecClusterVar(c)        (c)->var

/*! \brief create a vector cluster */
__ZM_EXPORT zVecCluster *zVecClusterCreate(zVecCluster *c, int core_size);

/*! \brief destroy a vector cluster */
__ZM_EXPORT void zVecClusterDestroy(zVecCluster *c);

/*! \brief the maximum silhouette in a cluster. */
__ZM_EXPORT double zVecClusterMaxSilhouette(const zVecCluster *c);

/*! \brief print a vector cluster to a file */
__ZM_EXPORT void zVecClusterFPrint(FILE *fp, const zVecCluster *c);
__ZM_EXPORT void zVecClusterValueFPrint(FILE *fp, const zVecCluster *c);

/* ********************************************************** */
/*! \brief methods for core and error computation
 *//* ******************************************************* */

ZDEF_STRUCT( __ZM_CLASS_EXPORT, zVecClusterMethod ){
  /*! \cond */
  zVec error; /* error vector */
  zVec (* error_fp)(const zVecClusterMethod*,const zVec,const zVec,void*,zVec); /* error function */
  void *error_util;
  double (* dist_fp)(const zVecClusterMethod*,const zVec,const zVec,void*); /* distance function */
  void *dist_util;
  int core_size; /* size of the core vector */
  zVec (* core_fp)(const zVecClusterMethod*,const zVecAddrList*,void*,zVec); /* a function to find core */
  void *core_util;
  /* for GMM */
  zVec (* lm_fp)(const zVecClusterMethod*,const zVecAddrList*,const double[],double,void*,zVec); /* loaded mean function */
  void *lm_util;
  zVec _lerr; /* workspace for loaded error vector */
  /*! \endcond */
};

#define zVecClusterMethodErrorF(cm,v1,v2) (cm)->error_fp( cm, (v1), (v2), (cm)->error_util, (cm)->error )
#define zVecClusterMethodDistF(cm,v1,v2)  (cm)->dist_fp( cm, (v1), (v2), (cm)->dist_util )
#define zVecClusterMethodCoreF(cm,vl,c)   (cm)->core_fp( cm, (vl), (cm)->core_util, (c) )
#define zVecClusterMethodLoadedMeanF(cm,vl,l,nk,m) (cm)->lm_fp( cm, (vl), (l), (nk), (cm)->lm_util, (m) )

/* methods for clustering based on linear-sum. */
__ZM_EXPORT zVec zVecClusterErrorLS(const zVecClusterMethod *method, const zVec p, const zVec core, void *dummy, zVec err);
__ZM_EXPORT zVec zVecClusterCoreLS(const zVecClusterMethod *method, const zVecAddrList *pl, void *dummy, zVec core);
__ZM_EXPORT zVec zVecClusterLoadedMeanLS(const zVecClusterMethod *method, const zVecAddrList *pl, const double load[], double n, void *dummy, zVec mean);

/*! \brief initialize methods for clustering. */
__ZM_EXPORT void zVecClusterMethodInit(zVecClusterMethod *method);

/*! \brief set default methods for clustering. */
__ZM_EXPORT zVecClusterMethod *zVecClusterMethodSetDefault(zVecClusterMethod *method, int size);

/*! \brief set a function to find core of a cluster. */
__ZM_EXPORT zVecClusterMethod *zVecClusterMethodSetCoreFunc(zVecClusterMethod *method, int core_size, zVec (* fp)(const zVecClusterMethod*,const zVecAddrList*,void*,zVec), void *util);
/*! \brief set a function to compute error of a sample in a cluster. */
__ZM_EXPORT zVecClusterMethod *zVecClusterMethodSetErrorFunc(zVecClusterMethod *method, int error_size, zVec (* fp)(const zVecClusterMethod*,const zVec,const zVec,void*,zVec), void *util);
/*! \brief set a function to compute the distance between two samples in a cluster. */
__ZM_EXPORT zVecClusterMethod *zVecClusterMethodSetDistFunc(zVecClusterMethod *method, double (* fp)(const zVecClusterMethod*,const zVec,const zVec,void*), void *util);
/*! \brief set a function to compute the loaded mean of samples (for GMM). */
__ZM_EXPORT zVecClusterMethod *zVecClusterMethodSetLoadedMeanFunc(zVecClusterMethod *method, zVec (* mean_l_fp)(const zVecClusterMethod*,const zVecAddrList*,const double[],double,void*,zVec), void *util);

/*! \brief set methods of linear-sum for clustering. */
__ZM_EXPORT zVecClusterMethod *zVecClusterMethodSetLS(zVecClusterMethod *method, void *util);

/*! \brief copy methods for clustering. */
__ZM_EXPORT zVecClusterMethod *zVecClusterMethodCopy(const zVecClusterMethod *src, zVecClusterMethod *dest);

/*! \brief destroy methods for clustering. */
__ZM_EXPORT void zVecClusterMethodDestroy(zVecClusterMethod *method);

/* ********************************************************** */
/*! \brief multiple vecter clusters class.
 *//* ******************************************************* */

zListClass( zVecClusterList, zVecClusterListCell, zVecCluster );

ZDEF_STRUCT( __ZM_CLASS_EXPORT, zVecMCluster ){
  zVecClusterList clist;    /*!< a list of clusters */
  zVecClusterMethod method; /*!< methods for core and error computation */
};

#define zVecMClusterClusterList(mc) ( &(mc)->clist )

/*! \brief initialize multiple vector clusters.
 *
 * zVecMClusterInit() intializes a multiple clusters \a mc.
 * \a size is the size of core vector, i.e., adjoint vector to samples, of a cluster, which is also
 * assigned to the size of error vector.
 * \return
 * zVecMClusterInit() returns a pointer \a mc if it succeeds to initialize \a mc.
 * Otherwise, the null pointer is returned.
 */
__ZM_EXPORT zVecMCluster *zVecMClusterInit(zVecMCluster *mc, int size);

#define zVecMClusterErrorF(mc,v1,v2) zVecClusterMethodErrorF( &(mc)->method, v1, v2 )
#define zVecMClusterDistF(mc,v1,v2)  zVecClusterMethodDistF( &(mc)->method, v1, v2 )
#define zVecMClusterCoreF(mc,vl,c)   zVecClusterMethodCoreF( &(mc)->method, vl, c )

#define zVecMClusterSetErrorFunc(mc,se,fp,util) zVecClusterMethodSetErrorFunc( &(mc)->method, (se), (fp), (util) )
#define zVecMClusterSetDistFunc(mc,fp,util)     zVecClusterMethodSetDistFunc( &(mc)->method, (fp), (util) )
#define zVecMClusterSetCoreFunc(mc,sc,fp,util)  zVecClusterMethodSetCoreFunc( &(mc)->method, (sc), (fp), (util) )

#define zVecMClusterSetLS(mc,util)   zVecClusterMethodSetLS( &(mc)->method, (util) )

/*! \brief allocate a set of vector clusters.
 *
 * zVecMClusterAlloc() allocates memory for multiple clusters \a mc.
 * \a n is the number of clusters.
 * \return
 * zVecMClusterAlloc() returns the pointer \a mc if it succeeds to allocate
 * memory. Otherwise, the null pointer is returned.
 */
__ZM_EXPORT zVecMCluster *zVecMClusterAlloc(zVecMCluster *mc, int n);

/*! \brief copy methods of a multiple cluster to another. */
#define zVecMClusterMethodCopy(src,dest) zVecClusterMethodCopy( &(src)->method, &(dest)->method )

/*! \brief move a multiple cluster to another. */
__ZM_EXPORT bool zVecMClusterMove(zVecMCluster *src, zVecMCluster *dest);

/*! \brief destroy a set of vector clusters. */
__ZM_EXPORT void zVecMClusterDestroy(zVecMCluster *mc);

/*! \brief evenness of clusters.
 *
 * zVecMClusterEvenness() returns the measure of evenness, namely, the ratio of
 * the maximum number of clustered samples over the minimum, of a multiple
 * cluster \a mc. \a mc is more even if it has as close value to 1.
 * \return
 * zVecMClusterEvenness() returns the measure of evenness.
 */
__ZM_EXPORT double zVecMClusterEvenness(const zVecMCluster *mc);

/*! \brief print vectors in a set of vector clusters */
__ZM_EXPORT void zVecMClusterFPrint(FILE *fp, const zVecMCluster *mc);
__ZM_EXPORT void zVecMClusterValueFPrint(FILE *fp[], const zVecMCluster *mc);
__ZM_EXPORT void zVecMClusterCoreFPrint(FILE *fp[], const zVecMCluster *mc);

/*! \brief print vectors in a set of clusters to files with a common basename. */
__ZM_EXPORT bool zVecMClusterValuePrintFile(const zVecMCluster *mc, const char *basename);

/* ********************************************************** */
/*! \brief clustering based on K-means
 *//* ******************************************************* */

/*! \brief clustering of vectors by K-means.
 *
 * zVecMClusterKMeansKKZ() and zVecMClusterKMeans() make multiple clusters from a set of
 * samples based on K-means method. zVecMClusterKMeansKKZ() uses KKZ method for
 * initialization, while zVecMClusterKMeans() uses K-means++.
 * For both functions, the resulted clusters are stored in \a mc. \a points is the
 * list of pointers to the original samples to be clusters. \a k is the number of
 * clusters.
 * \return
 * zVecMClusterKMeansKKZ() and zVecMClusterKMeans() return the number of iterations that
 * is taken in the K-means method.
 */
__ZM_EXPORT int zVecMClusterKMeansKKZ(zVecMCluster *mc, const zVecAddrList *points, int k);
__ZM_EXPORT int zVecMClusterKMeans(zVecMCluster *mc, const zVecAddrList *points, int k);

/*! \brief clustering of vectors by K-medoids.
 *
 * zVecMClusterKMedoids() makes multiple clusters from a set of samples based on K-medoids
 * method. It initializes clusters by the same method with K-means++.
 * The resulted clusters are stored in \a mc. \a points is the list of pointers to
 * the original samples to be clusters. \a k is the number of clusters.
 * \return
 * zVecMClusterKMedoids() returns the number of iterations that is taken in the method.
 * \sa
 * zVecMClusterKMeansKKZ, zVecMClusterKMeans
 */
__ZM_EXPORT int zVecMClusterKMedoids(zVecMCluster *mc, const zVecAddrList *points, int k);

/*! \brief compute the mean silhouette of a set of clusters. */
__ZM_EXPORT double zVecMClusterMeanSilhouette(zVecMCluster *mc);

/*! \brief print silhouettes of a set of vector clusters to files. */
__ZM_EXPORT bool zVecMClusterSilhouettePrintFile(const zVecMCluster *mc, const char *basename);

/* ********************************************************** */
/*! \brief clustering based on X-means
 *//* ******************************************************* */

/*! \brief clustering of vectors by X-means based on hyperdensity of clusters. */
__ZM_EXPORT int zVecMClusterXMeansDensity(zVecMCluster *mc, const zVecAddrList *points);

/*! \brief clustering of vectors by X-means based on BIC */
__ZM_EXPORT int zVecMClusterXMeansBIC(zVecMCluster *mc, const zVecAddrList *points);

__END_DECLS

#endif /* __ZM_MVA_CLUSTER_H__ */
