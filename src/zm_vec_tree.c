/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_tree - vector binary tree (kd-tree) class.
 */

#include <zm/zm_vec.h>

/* initialize a 3D vector tree. */
zVecTree *zVecTreeInit(zVecTree *tree, int dim)
{
  tree->split = -1; /* invalid split index */
  tree->s[0] = tree->s[1] = NULL;
  tree->v = NULL;
  tree->vmin = zVecAlloc( dim );
  tree->vmax = zVecAlloc( dim );
  if( !tree->vmin || !tree->vmax ){
    zVecFree( tree->vmin );
    zVecFree( tree->vmax );
    return NULL;
  }
  zVecSetAll( tree->vmin,-HUGE_VAL );
  zVecSetAll( tree->vmax, HUGE_VAL );
  return tree;
}

/* destroy a 3D vector tree. */
void zVecTreeDestroy(zVecTree *tree)
{
  if( tree->s[0] ){
    zVecTreeDestroy( tree->s[0] );
    free( tree->s[0] );
  }
  if( tree->s[1] ){
    zVecTreeDestroy( tree->s[1] );
    free( tree->s[1] );
  }
  zVecFree( tree->v );
  zVecFree( tree->vmin );
  zVecFree( tree->vmax );
}

/* create a leaf of a 3D vector tree. */
static zVecTree *_zVecTreeCreateLeaf(int split, zVec v)
{
  zVecTree *leaf;

  if( !( leaf = zAlloc( zVecTree, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !( leaf->v = zVecClone( v ) ) ){
    ZALLOCERROR();
    free( leaf );
    return NULL;
  }
  leaf->split = split;
  leaf->s[0] = leaf->s[1] = NULL;
  return leaf;
}

/* return an index of a node which contains a given 3D vector. */
static int _zVecTreeChooseBranch(const zVecTree *node, zVec v)
{
  return zVecElemNC(v,node->split) >= zVecElemNC(node->v,node->split) ? 0 : 1;
}

/* add a new 3D vector to a tree. */
static zVecTree *_zVecTreeAdd(zVecTree *node, zVec v)
{
  int b;
  zVecTree *leaf;

  if( node->s[( b = _zVecTreeChooseBranch( node, v ) )] )
    return _zVecTreeAdd( node->s[b], v );
  if( !( leaf = _zVecTreeCreateLeaf( ( node->split + 1 ) % zVecSizeNC(v), v ) ) )
    return NULL;
  node->s[b] = leaf;
  leaf->vmin = zVecClone( node->vmin );
  leaf->vmax = zVecClone( node->vmax );
  if( !leaf->vmin || !leaf->vmax ){
    zVecFree( leaf->vmin );
    zVecFree( leaf->vmax );
    free( leaf );
    return NULL;
  }
  if( b == 0 )
    zVecElemNC(leaf->vmin,node->split) = zVecElemNC(node->v,node->split);
  else /* b == 1 */
    zVecElemNC(leaf->vmax,node->split) = zVecElemNC(node->v,node->split);
  return leaf;
}

/* add a new 3D vector to a tree. */
zVecTree *zVecTreeAdd(zVecTree *tree, zVec v)
{
  if( tree->split == -1 ){
    tree->split = 0;
    return ( tree->v = zVecClone( v ) ) ? tree : NULL;
  }
  return _zVecTreeAdd( tree, v );
}

/* nearest neighbor search */

/* check if a sphere is overlapped with a bounding box of a node. */
static bool _zVecTreeIsOverlap(const zVecTree *node, const zVec c, double r)
{
  int i;
  double d, vd;

  for( d=0, i=0; i<zVecSizeNC(c); i++ ){
    vd = zVecElemNC(c,i);
    if( vd < zVecElemNC(node->vmin,i) )
      d += zSqr( vd - zVecElemNC(node->vmin,i) );
    if( vd > zVecElemNC(node->vmax,i) )
      d += zSqr( vd - zVecElemNC(node->vmax,i) );
  }
  return d <= r*r + zTOL ? true : false;
}

/* test if a node is the current nearest neighbor to a 3D vector. */
static void _zVecTreeNNTest(const zVecTree *node, const zVec v, zVecTree **nn, double *dmin)
{
  double d;

  if( ( d = zVecDist( node->v, v ) ) < *dmin ){
    *nn = (zVecTree *)node;
    *dmin = d;
  }
}

/* an internal recursive call of the nearest neighbor search;
 * check the opposite side of branch.
 */
static double _zVecTreeNNOpp(const zVecTree *node, const zVec v, zVecTree **nn, double *dmin)
{
  _zVecTreeNNTest( node, v, nn, dmin );
  if( node->s[0] && _zVecTreeIsOverlap( node->s[0], v, *dmin ) )
    _zVecTreeNNOpp( node->s[0], v, nn, dmin );
  if( node->s[1] && _zVecTreeIsOverlap( node->s[1], v, *dmin ) )
    _zVecTreeNNOpp( node->s[1], v, nn, dmin );
  return *dmin;
}

/* an internal recursive call of the nearest neighbor search. */
static double _zVecTreeNN(const zVecTree *node, const zVec v, zVecTree **nn, double *dmin)
{
  int b;
  zVecTree *ob;

  if( node->s[( b = _zVecTreeChooseBranch( node, v ) )] )
    _zVecTreeNN( node->s[b], v, nn, dmin );
  _zVecTreeNNTest( node, v, nn, dmin );
  if( ( ob = node->s[1-b] ) && _zVecTreeIsOverlap( ob, v, *dmin ) )
    _zVecTreeNNOpp( ob, v, nn, dmin );
  return *dmin;
}

/* find the nearest neighbor to a 3D vector in a tree. */
double zVecTreeNN(const zVecTree *tree, const zVec v, zVecTree **nn)
{
  double dmin = HUGE_VAL;

  *nn = NULL;
  return _zVecTreeNN( tree, v, nn, &dmin );
}
