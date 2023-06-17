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

/* find the nearest neighbor of a vector by a naive algorithm. */
zVec zVecListNN(zVecList *list, zVec v, double *dmin)
{
  zVecListCell *cell;
  double d, __dmin;
  zVec nn = NULL;

  if( dmin == NULL ) dmin = &__dmin; /* dummy pointer */
  *dmin = HUGE_VAL;
  zListForEach( list, cell )
    if( ( d = zVecSqrDist( cell->data, v ) ) < *dmin ){
      *dmin = d;
      nn = cell->data;
    }
  *dmin = sqrt( *dmin );
  return nn;
}

/* scan vectors from a file and creates a list of them. */
zVecList *zVecListFScan(FILE *fp, zVecList *list)
{
  zVec v;

  zListInit( list );
  while( ( v = zVecDataFScan( fp ) ) )
    zVecListInsertHead( list, v );
  if( zListIsEmpty( list ) )
    ZRUNWARN( ZM_WARN_VECLIST_EMPTY );
  return list;
}

/* print vectors in a list to a file. */
void zVecListFPrint(FILE *fp, zVecList *list)
{
  zVecListCell *cp;

  zListForEach( list, cp )
    zVecDataFPrint( fp, cp->data );
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
