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

/*! \struct zVecCluster
 * \brief vecter cluster class.
 */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zVecCluster ){
  zVecAddrList vlist; /*!< \brief list of vectors */
  zVec core;          /*!< \brief core vector of the cluster */
  double var;         /*!< \brief variance */
  /*! \cond */
  double *_sil;       /* silhouettes */
  /*! \endcond */
};

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

/*! \struct zVecClusterMethod
 * \brief methods for core and error computation.
 */
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

/*! \struct zVecClusterList
 * \brief list of vecter clusters.
 */
ZEDA_DEF_LIST_CLASS( zVecClusterList, zVecClusterListCell, zVecCluster );

/*! \struct zVecMultiCluster
 * \brief multiple vector clusters.
 */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zVecMultiCluster ){
  zVecClusterList clist;    /*!< a list of clusters */
  zVecClusterMethod method; /*!< methods for core and error computation */
};

#define zVecMultiClusterClusterList(mc) ( &(mc)->clist )

/*! \brief initialize multiple vector clusters.
 *
 * zVecMultiClusterInit() intializes a multiple clusters \a mc.
 * \a size is the size of core vector, i.e., adjoint vector to samples, of a cluster, which is also
 * assigned to the size of error vector.
 * \return
 * zVecMultiClusterInit() returns a pointer \a mc if it succeeds to initialize \a mc.
 * Otherwise, the null pointer is returned.
 */
__ZM_EXPORT zVecMultiCluster *zVecMultiClusterInit(zVecMultiCluster *mc, int size);

#define zVecMultiClusterErrorF(mc,v1,v2) zVecClusterMethodErrorF( &(mc)->method, v1, v2 )
#define zVecMultiClusterDistF(mc,v1,v2)  zVecClusterMethodDistF( &(mc)->method, v1, v2 )
#define zVecMultiClusterCoreF(mc,vl,c)   zVecClusterMethodCoreF( &(mc)->method, vl, c )

#define zVecMultiClusterSetErrorFunc(mc,se,fp,util) zVecClusterMethodSetErrorFunc( &(mc)->method, (se), (fp), (util) )
#define zVecMultiClusterSetDistFunc(mc,fp,util)     zVecClusterMethodSetDistFunc( &(mc)->method, (fp), (util) )
#define zVecMultiClusterSetCoreFunc(mc,sc,fp,util)  zVecClusterMethodSetCoreFunc( &(mc)->method, (sc), (fp), (util) )

#define zVecMultiClusterSetLS(mc,util)   zVecClusterMethodSetLS( &(mc)->method, (util) )

/*! \brief allocate a set of vector clusters.
 *
 * zVecMultiClusterAlloc() allocates memory for multiple clusters \a mc.
 * \a n is the number of clusters.
 * \return
 * zVecMultiClusterAlloc() returns the pointer \a mc if it succeeds to allocate
 * memory. Otherwise, the null pointer is returned.
 */
__ZM_EXPORT zVecMultiCluster *zVecMultiClusterAlloc(zVecMultiCluster *mc, int n);

/*! \brief copy methods of a multiple cluster to another. */
#define zVecMultiClusterMethodCopy(src,dest) zVecClusterMethodCopy( &(src)->method, &(dest)->method )

/*! \brief move a multiple cluster to another. */
__ZM_EXPORT bool zVecMultiClusterMove(zVecMultiCluster *src, zVecMultiCluster *dest);

/*! \brief destroy a set of vector clusters. */
__ZM_EXPORT void zVecMultiClusterDestroy(zVecMultiCluster *mc);

/*! \brief evenness of clusters.
 *
 * zVecMultiClusterEvenness() returns the measure of evenness, namely, the ratio of
 * the maximum number of clustered samples over the minimum, of a multiple
 * cluster \a mc. \a mc is more even if it has as close value to 1.
 * \return
 * zVecMultiClusterEvenness() returns the measure of evenness.
 */
__ZM_EXPORT double zVecMultiClusterEvenness(const zVecMultiCluster *mc);

/*! \brief print vectors in a set of vector clusters */
__ZM_EXPORT void zVecMultiClusterFPrint(FILE *fp, const zVecMultiCluster *mc);
__ZM_EXPORT void zVecMultiClusterValueFPrint(FILE *fp[], const zVecMultiCluster *mc);
__ZM_EXPORT void zVecMultiClusterCoreFPrint(FILE *fp[], const zVecMultiCluster *mc);

/*! \brief print vectors in a set of clusters to files with a common basename. */
__ZM_EXPORT bool zVecMultiClusterValuePrintFile(const zVecMultiCluster *mc, const char *basename);

/* ********************************************************** */
/*! \brief clustering based on K-means
 *//* ******************************************************* */

/*! \brief clustering of vectors by K-means.
 *
 * zVecMultiClusterKMeansKKZ() and zVecMultiClusterKMeans() make multiple clusters from a set of
 * samples based on K-means method. zVecMultiClusterKMeansKKZ() uses KKZ method for
 * initialization, while zVecMultiClusterKMeans() uses K-means++.
 * For both functions, the resulted clusters are stored in \a mc. \a points is the
 * list of pointers to the original samples to be clusters. \a k is the number of
 * clusters.
 * \return
 * zVecMultiClusterKMeansKKZ() and zVecMultiClusterKMeans() return the number of iterations that
 * is taken in the K-means method.
 */
__ZM_EXPORT int zVecMultiClusterKMeansKKZ(zVecMultiCluster *mc, const zVecAddrList *points, int k);
__ZM_EXPORT int zVecMultiClusterKMeans(zVecMultiCluster *mc, const zVecAddrList *points, int k);

/*! \brief clustering of vectors by K-medoids.
 *
 * zVecMultiClusterKMedoids() makes multiple clusters from a set of samples based on K-medoids
 * method. It initializes clusters by the same method with K-means++.
 * The resulted clusters are stored in \a mc. \a points is the list of pointers to
 * the original samples to be clusters. \a k is the number of clusters.
 * \return
 * zVecMultiClusterKMedoids() returns the number of iterations that is taken in the method.
 * \sa
 * zVecMultiClusterKMeansKKZ, zVecMultiClusterKMeans
 */
__ZM_EXPORT int zVecMultiClusterKMedoids(zVecMultiCluster *mc, const zVecAddrList *points, int k);

/*! \brief compute the mean silhouette of a set of clusters. */
__ZM_EXPORT double zVecMultiClusterMeanSilhouette(zVecMultiCluster *mc);

/*! \brief print silhouettes of a set of vector clusters to files. */
__ZM_EXPORT bool zVecMultiClusterSilhouettePrintFile(const zVecMultiCluster *mc, const char *basename);

/* ********************************************************** */
/*! \brief clustering based on X-means
 *//* ******************************************************* */

/*! \brief clustering of vectors by X-means based on hyperdensity of clusters. */
__ZM_EXPORT int zVecMultiClusterXMeansDensity(zVecMultiCluster *mc, const zVecAddrList *points);

/*! \brief clustering of vectors by X-means based on BIC */
__ZM_EXPORT int zVecMultiClusterXMeansBIC(zVecMultiCluster *mc, const zVecAddrList *points);

__END_DECLS

#endif /* __ZM_MVA_CLUSTER_H__ */
