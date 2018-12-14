/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_rrt.h
 * \brief RRT
 * \author Zhidao
 */

#ifndef __ZM_RRT_H__
#define __ZM_RRT_H__

#include <zm/zm_vec.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup RRT
 * \{ *//* ************************************************** */

/*! \brief RRT node class to construct open tree
 */
typedef struct _zRRTNode{
  struct _zRRTNode *parent; /*!< a pointer to the parent */
  zVec v;                   /*! data vector */
} zRRTNode;

/*! \brief RRT node list */
zListClass( zRRTList, zRRTListCell, zRRTNode );

/*! \brief RRT solver */
typedef struct{
  zRRTList slist; /*!< RRT from start node */
  zRRTList glist; /*!< RRT from goal node */
  zVec min;       /*!< minimum bound of vectors */
  zVec max;       /*!< maximum bound of vectors */
  double eps;     /*!< delta length for extend-tree operation */

  double (* dist)(zVec,zVec,void*); /* distance function to define the metric */
  zVec (* ext)(zVec,zVec,double,zVec,void*); /* extension function */
  bool (* chk)(zVec,void*); /* feasibility checking function */
} zRRT;

__EXPORT void zRRTInit(zRRT *rrt, zVec min, zVec max, double eps, double (* dist)(zVec,zVec,void*), zVec (* ext)(zVec,zVec,double,zVec,void*), bool (* chk)(zVec,void*));
__EXPORT void zRRTDestroy(zRRT *rrt);

__EXPORT bool zRRTConnect(zRRT *rrt, zVec start, zVec goal, int iter, void *util, zVecList *path);
__EXPORT void zRRTPathShortcut(zRRT *rrt, void *util, zVecList *path);

#define zRRTPathDestroy(path) zVecListDestroy( path, true )

/*! \brief RRT-Escapement solver (proposed by Y. Shimizu in 2012) */
__EXPORT bool zRRTEsc(zRRT *rrt, zVec start, int iter, void *util, zVec goal);

/*! \} */

__END_DECLS

#endif /* __ZM_RRT_H__ */
