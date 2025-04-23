/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_list - vector list class.
 */

#include <zm/zm_vec.h>

/* create a vector list cell. */
static zVecListCell *_zVecListCellCreate(zVec v, bool flag)
{
  zVecListCell *cell;

  if( !( cell = zAlloc( zVecListCell, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( flag ){
    if( !( cell->data = zVecClone( v ) ) ){
      zFree( cell );
      return NULL;
    }
  } else
    cell->data = v;
  return cell;
}

/* insert a vector list cell at the head of a list. */
static zVecListCell *_zVecListInsertHead(zVecList *list, zVec v, bool flag)
{
  zVecListCell *cell;

  if( ( cell = _zVecListCellCreate( v, flag ) ) )
    zListInsertHead( list, cell );
  return cell;
}

/* insert a vector list cell at the tail of a list. */
static zVecListCell *_zVecListInsertTail(zVecList *list, zVec v, bool flag)
{
  zVecListCell *cell;

  if( ( cell = _zVecListCellCreate( v, flag ) ) )
    zListInsertTail( list, cell );
  return cell;
}

/* ********************************************************** */
/* CLASS: zVecList
 * vector list class.
 * ********************************************************** */

/* insert a vector list cell at the head of a list. */
zVecListCell *zVecListInsertHead(zVecList *list, zVec v)
{
  return _zVecListInsertHead( list, v, true );
}

/* insert a vector list cell at the tail of a list. */
zVecListCell *zVecListInsertTail(zVecList *list, zVec v)
{
  return _zVecListInsertTail( list, v, true );
}

/* destroy a vector list. */
void zVecListDestroy(zVecList *list)
{
  zVecListCell *cell;

  while( !zListIsEmpty( list ) ){
    zListDeleteHead( list, &cell );
    zFree( cell->data );
    zFree( cell );
  }
}

/* select a vector of a vector list randomly. */
zVec zVecListSelectRand(zVecList *vl)
{
  zVecListCell *cp;
  int i;

  i = zRandI( 0, zListSize(vl)-1 );
  zListItem( vl, i, &cp );
  return cp->data;
}

/* find the nearest neighbor of a vector by a naive algorithm. */
double zVecListNN(const zVecList *list, const zVec v, zVec *nn)
{
  zVecListCell *cell;
  double d2, dmin2;

  *nn = NULL;
  dmin2 = HUGE_VAL;
  zListForEach( list, cell )
    if( ( d2 = zVecSqrDist( cell->data, v ) ) < dmin2 ){
      dmin2 = d2;
      *nn = cell->data;
    }
  return sqrt( dmin2 );
}

/* default distance function for Ramer-Douglas-Peucker algorithm. */
static double _zVecEdgeDistRDP(const zVec v, const zVec v1, const zVec v2, void *dummy)
{
  return zVecEdgeDist( v, v1, v2 );
}

/* internal recursive function or Ramer-Douglas-Peucker algorithm. */
static void _zVecListRDP(zVecList *list, zVecListCell *tail, zVecListCell *head, double (* metric)(const zVec,const zVec,const zVec,void*), void *util, double tol)
{
  zVecListCell *cp, *cp_farthest = NULL;
  double d, d_max = 0;

  for( cp=zListCellNext(tail); cp!=head; cp=zListCellNext(cp) ){
    if( ( d = metric( cp->data, tail->data, head->data, util ) ) > d_max ){
      d_max = d;
      cp_farthest = cp;
    }
  }
  if( d_max > tol ){
    _zVecListRDP( list, tail, cp_farthest, metric, util, tol );
    _zVecListRDP( list, cp_farthest, head, metric, util, tol );
  } else{
    while( ( cp = zListCellNext(tail) ) != head ){
      zListPurge( list, cp );
      zVecFree( cp->data );
      free( cp );
    }
  }
}

/* Ramer-Douglas-Peucker algorithm. */
void zVecListRDP(zVecList *list, double (* metric)(const zVec,const zVec,const zVec,void*), void *util, double tol)
{
  if( !metric ) metric = _zVecEdgeDistRDP;
  _zVecListRDP( list, zListTail(list), zListHead(list), metric, util, tol );
}

/* scan vectors from a file and creates a list of them. */
zVecList *zVecListFScan(FILE *fp, zVecList *list)
{
  zVec v;

  zListInit( list );
  while( ( v = zVecValueFScan( fp ) ) )
    zVecListInsertHead( list, v );
  if( zListIsEmpty( list ) )
    ZRUNWARN( ZM_WARN_VECLIST_EMPTY );
  return list;
}

/* print vectors in a list to a file. */
void zVecListFPrint(FILE *fp, const zVecList *list)
{
  zVecListCell *cp;

  zListForEach( list, cp )
    zVecValueFPrint( fp, cp->data );
}

/* ********************************************************** */
/* CLASS: zVecAddrList
 * vector address list class.
 * ********************************************************** */

/* insert a pointer to a vector to a list at the head. */
zVecListCell *zVecAddrListInsertHead(zVecAddrList *list, zVec v)
{
  return _zVecListInsertHead( list, v, false );
}

/* insert a pointer to a vector to a list at the tail. */
zVecListCell *zVecAddrListInsertTail(zVecAddrList *list, zVec v)
{
  return _zVecListInsertTail( list, v, false );
}

/* destroy a vector list. */
void zVecAddrListDestroy(zVecList *list)
{
  zVecAddrListCell *cell;

  while( !zListIsEmpty( list ) ){
    zListDeleteHead( list, &cell );
    zFree( cell );
  }
}
