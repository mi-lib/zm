/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_tree - vector binary tree (kd-tree) class.
 */

#ifndef __ZM_VEC_TREE_H__
#define __ZM_VEC_TREE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zVecTree
 * vector binary tree class.
 * ********************************************************** */

typedef struct _zVecTree{
  int split; /*!< split index */
  zVec v;    /*!< spliting vertex */
  zVec vmin; /*!< minimum corner of bounding box */
  zVec vmax; /*!< maximum corner of bounding box */
  struct _zVecTree *s[2]; /*!< binary branches */
} zVecTree;

__EXPORT zVecTree *zVecTreeInit(zVecTree *tree, int dim);
__EXPORT void zVecTreeDestroy(zVecTree *tree);

__EXPORT zVecTree *zVecTreeAdd(zVecTree *tree, zVec v);

__EXPORT double zVecTreeNN(zVecTree *tree, zVec v, zVecTree **nn);

__END_DECLS

#endif /* __ZM_VEC_TREE_H__ */
