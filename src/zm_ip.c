/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip - interpolation.
 */

#include <zm/zm_ip.h>

/* ********************************************************** */
/* CLASS: zIP
 * interpolation curve class
 * ********************************************************** */

/* zIPDataAlloc
 * - assign a sequence to an interpolator and allocate
 *   the internal workspace.
 */
bool zIPDataAlloc(zIPData *dat, zSeq *seq)
{
  zSeqListCell *cp;
  zIPKnotCell *kp;
  double t = 0;
  register int i = 0;

  if( zListIsEmpty(seq) ){
    ZRUNWARN( ZM_WARN_SEQ_EMPTY );
    return true;
  }
  zArrayAlloc( &dat->knot, zIPKnotCell, zListNum(seq) );
  if( zArrayNum(&dat->knot) == 0 ) return false;
  if( !zVecArrayAlloc( &dat->va, zVecSize(zListHead(seq)->data.v), zListNum(seq) ) ){
    ZALLOCERROR();
    zIPDataFree( dat );
    return false;
  }
  dat->seq = seq;
  zListForEachRew( dat->seq, cp ){
    kp = zIPKnot(dat,i++);
    kp->t = ( t += cp->data.dt );
    kp->cp = cp;
  }
  return true;
}

/* zIPDataFree
 * - free the internal workspace of an interpolator.
 */
void zIPDataFree(zIPData *dat)
{
  zArrayFree( &dat->knot );
  zVecArrayFree( &dat->va );
  dat->seq = NULL;
}

/* zIPSeg
 * - find a segment in which the absessi is included.
 */
int zIPSeg(zIPData *dat, double t)
{
  register int i, j, k;

  for( i=0, j=zIPSize(dat)-1; ; ){
    if( ( k = ( i + j ) / 2 ) == i ) break;
    if( zIPTime(dat,k) > t ) j = k;
    else i = k;
  }
  return i;
}
