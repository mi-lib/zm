/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_tree - vector binary tree (kd-tree) class.
 */

#ifndef __ZM_VEC_TREE_H__
#define __ZM_VEC_TREE_H__

/* NOTE: never include this header file in user programs. */

/*! \struct zVecTree
 * \brief binary tree of vectors.
 */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zVecTree ){
  int split; /*!< split index */
  zVec v;    /*!< spliting vertex */
  zVec vmin; /*!< minimum corner of bounding box */
  zVec vmax; /*!< maximum corner of bounding box */
  zVecTree *s[2]; /*!< binary branches */
};

__BEGIN_DECLS

__ZM_EXPORT zVecTree *zVecTreeInit(zVecTree *tree, int dim);
__ZM_EXPORT void zVecTreeDestroy(zVecTree *tree);

__ZM_EXPORT zVecTree *zVecTreeAdd(zVecTree *tree, zVec v);

__ZM_EXPORT double zVecTreeNN(const zVecTree *tree, const zVec v, zVecTree **nn);

__END_DECLS

#endif /* __ZM_VEC_TREE_H__ */
