/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_list - vector list class.
 */

#include <zm/zm_vec.h>

static zVecListCell *_zVecListCellCreate(zVec v, bool flag);

/* (static)
 * _zVecListCellCreate
 * - create a vector list cell.
 */
zVecListCell *_zVecListCellCreate(zVec v, bool flag)
{
  zVecListCell *cell;

  if( !( cell = zAlloc( zVecListCell, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( flag ){
    if( !( cell->data = zVecClone( v ) ) ){
      ZALLOCERROR();
      zFree( cell );
      return NULL;
    }
  } else{
    cell->data = v;
  }
  return cell;
}

/* zVecListInsertHead
 * - insert a vector list cell at the head of a list.
 */
zVecListCell *zVecListInsertHead(zVecList *list, zVec v, bool flag)
{
  zVecListCell *cell;

  if( ( cell = _zVecListCellCreate( v, flag ) ) )
    zListInsertHead( list, cell );
  return cell;
}

/* zVecListInsertTail
 * - insert a vector list cell at the tail of a list.
 */
zVecListCell *zVecListInsertTail(zVecList *list, zVec v, bool flag)
{
  zVecListCell *cell;

  if( ( cell = _zVecListCellCreate( v, flag ) ) )
    zListInsertTail( list, cell );
  return cell;
}

/* zVecListDestroy
 * - destroy a vector list.
 */
void zVecListDestroy(zVecList *list, bool flag)
{
  zVecListCell *cell;

  if( flag ){
    while( !zListIsEmpty( list ) ){
      zListDeleteHead( list, &cell );
      zFree( cell->data );
      zFree( cell );
    }
  } else{
    while( !zListIsEmpty( list ) ){
      zListDeleteHead( list, &cell );
      zFree( cell );
    }
  }
}

/* zVecListNN
 * - find the nearest neighbor of a vector by a naive algorithm.
 */
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

/* zVecListFWrite
 * - output vectors in a list.
 */
void zVecListFWrite(FILE *fp, zVecList *list)
{
  zVecListCell *cp;

  zListForEach( list, cp )
    zVecDataFWrite( fp, cp->data );
}
