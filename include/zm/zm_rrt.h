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

/*! \brief RRT node class to construct open tree */
typedef struct _zRRTNode{
  zVec v;                   /*!< a point vector */
  double cost_edge;         /*!< edge cost to the parent */
  double cost_total;        /*!< total cost */
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
__ZM_EXPORT void zRRTInit(zRRT *rrt, zVec min, zVec max, double eps, double (* distance)(zVec,zVec,void*), zVec (* extend)(zVec,zVec,double,zVec,void*), bool (* chk_collision)(zVec,void*), bool (* chk_goal)(zVec,void*));
/*! \brief destroy RRTs. */
__ZM_EXPORT void zRRTDestroy(zRRT *rrt);

/*! \brief find a path based on the RRT algorithm.
 *
 * The original paper about the algorithm is:
 * Rapidly-exploring random trees: A new tool for path planning. S. M. LaValle. TR 98-11,
 * Computer Science Dept., Iowa State University, October 1998.
 */
__ZM_EXPORT bool zRRTFindPath(zRRT *rrt, zVec start, int iter, void *util, zVecList *path, double *cost);

/*! \brief find the optimum path based on the RRT* algorithm.
 *
 * The original paper about the algorithm is:
 * Sampling-based Algorithms for Optimal Motion Planning. Sertac Karaman and Emilio Frazzoli.
 * arXiv:1105.1186
 */
__ZM_EXPORT bool zRRTFindPathOpt(zRRT *rrt, zVec start, int iter, void *util, zVecList *path, double *cost);

/*! \brief find a path based on the RRT-connect algorithm.
 *
 * The original paper about the algorithm is:
 * RRT-connect: An efficient approach to single-query path planning. J. J. Kuffner and
 * S. M. LaValle. In Proceedings IEEE International Conference on Robotics and Automation,
 * pp. 995-1001, 2000.
 */
__ZM_EXPORT bool zRRTFindPathDual(zRRT *rrt, zVec start, zVec goal, int iter, void *util, zVecList *path, double *cost);

/*! \brief a postprocess for RRT family to shortcut a path. */
__ZM_EXPORT bool zRRTShortcutPath(zRRT *rrt, void *util, zVecList *path, double *cost);

/*! \brief RRT-escapement algorithm to find a collision-free point proposed by Y. Shimizu in 2012. */
__ZM_EXPORT bool zRRTEscape(zRRT *rrt, zVec start, int iter, void *util, zVec goal);

/*! \} */

__END_DECLS

#endif /* __ZM_RRT_H__ */
