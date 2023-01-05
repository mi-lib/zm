/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_rrt.h
 * \brief Rapidly-explored Random Tree algorithm and its family.
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
  zVec v;                   /*!< a point vector */
  struct _zRRTNode *parent; /*!< pointer to the parent node */
} zRRTNode;

/*! \brief RRT node list */
zListClass( zRRTList, zRRTListCell, zRRTNode );

/*! \brief data structure for RRT family */
typedef struct{
  zRRTList slist; /*!< RRT from start node */
  zRRTList glist; /*!< RRT from goal node (for RRT-connect) */
  zVec min;       /*!< minimum bound of vectors */
  zVec max;       /*!< maximum bound of vectors */
  double eps;     /*!< delta length for extend-tree operation */

  /*! \cond */
  double (* _distance)(zVec,zVec,void*); /* distance function to define metric */
  zVec (* _extend)(zVec,zVec,double,zVec,void*); /* extension function */
  bool (* _check_collision)(zVec,void*); /* collision-check function */
  bool (* _check_goal)(zVec,void*); /* goal-check function */
  /*! \endcond */
} zRRT;

/*! \brief set the distance function. */
#define zRRTSetDistanceFunc(rrt,f)       ( (rrt)->_distance = (f) )
/*! \brief set the extension function. */
#define zRRTSetExtendFunc(rrt,f)         ( (rrt)->_extend = (f) )
/*! \brief set the collision-check function. */
#define zRRTSetCheckCollisionFunc(rrt,f) ( (rrt)->_check_collision = (f) )
/*! \brief set the goal-check function (for the original RRT). */
#define zRRTSetCheckGoalFunc(rrt,f)      ( (rrt)->_check_goal = (f) )

/*! \brief initialize RRTs. */
__EXPORT void zRRTInit(zRRT *rrt, zVec min, zVec max, double eps, double (* distance)(zVec,zVec,void*), zVec (* extend)(zVec,zVec,double,zVec,void*), bool (* chk_collision)(zVec,void*), bool (* chk_goal)(zVec,void*));
/*! \brief destroy RRTs. */
__EXPORT void zRRTDestroy(zRRT *rrt);

/*! \brief find a path based on the RRT algorithm. */
__EXPORT bool zRRTFindPath(zRRT *rrt, zVec start, int iter, void *util, zVecList *path);
/*! \brief find a path based on the RRT-connect algorithm. */
__EXPORT bool zRRTFindPathDual(zRRT *rrt, zVec start, zVec goal, int iter, void *util, zVecList *path);

/*! \brief a postprocess for RRT family to shortcut a path. */
__EXPORT void zRRTShortcutPath(zRRT *rrt, void *util, zVecList *path);

/*! \brief RRT-escapement algorithm to find a collision-free point proposed by Y. Shimizu in 2012. */
__EXPORT bool zRRTEscape(zRRT *rrt, zVec start, int iter, void *util, zVec goal);

/*! \} */

__END_DECLS

#endif /* __ZM_RRT_H__ */
