/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_seq - vector sequence.
 */

#include <zm/zm_seq.h>

/* ********************************************************** */
/* CLASS: zSeq
 * motion sequence class for multiple dimension system.
 * ********************************************************** */

/* initialize a sequence. */
zSeq *zSeqInit(zSeq *seq)
{
  zListInit( seq );
  zListRoot(seq)->data.dt = 0;
  zListRoot(seq)->data.v = NULL;
  return seq;
}

/* free a sequence. */
void zSeqFree(zSeq *seq)
{
  zSeqCell *cp;

  while( !zListIsEmpty( seq ) )
    if( ( cp = zSeqDequeue( seq ) ) )
      zSeqCellFree( cp );
}

/* enqueue a vector to a sequence. */
zSeqCell *zSeqEnqueue(zSeq *seq, zVec v, double dt)
{
  zSeqCell *cp;

  if( !v ){
    ZRUNERROR( ZM_ERR_VEC_NULL );
    return NULL;
  }
  if( !( cp = zAlloc( zSeqCell, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  cp->data.v = v;
  cp->data.dt = dt;
  zQueueEnqueue( seq, cp );
  return cp;
}

/* dequeue a vector from a sequence. */
zSeqCell *zSeqDequeue(zSeq *seq)
{
  zSeqCell *cp;

  if( zListIsEmpty( seq ) ){
    ZRUNWARN( ZM_WARN_SEQ_EMPTY );
    return NULL;
  }
  zQueueDequeue( seq, &cp );
  return cp;
}

/* jump to the specified frame of a sequence. */
zSeqCell *zSeqJump(zSeq *seq, int step)
{
  zSeqCell *cp;

  for( cp=zListHead(seq); ; step--, cp=zListCellPrev(cp) ){
    if( cp == zListRoot(seq) ){
      ZRUNWARN( ZM_WARN_SEQ_STEP );
      return zListTail(seq);
    }
    if( step <= 0 ) break;
  }
  return cp;
}

/* scan a sequence from a file. */
bool zSeqScanFile(zSeq *seq, char filename[])
{
  FILE *fp;

  if( !( fp = zOpenFile( filename, ZSEQ_SUFFIX, "r" ) ) )
    return false;
  zSeqFScan( fp, seq );
  fclose( fp );
  return true;
}

/* scan a sequence from a file. */
zSeq *zSeqFScan(FILE *fp, zSeq *seq)
{
  zSeqCell *cp;

  zSeqInit( seq );
  while( !feof( fp ) ){
    if( !zFSkipDelimiter( fp ) ) break;
    if( !( cp = zAlloc( zSeqCell, 1 ) ) ){
      ZALLOCERROR();
      break;
    }
    if( !zFDouble( fp, &cp->data.dt ) ){
      ZRUNERROR( ZM_ERR_SEQ_DT_NOTFOUND );
      break;
    }
    if( !( cp->data.v = zVecFScan( fp ) ) ){
      ZALLOCERROR();
      free( cp );
      break;
    }
    zQueueEnqueue( seq, cp );
  }
  return seq;
}

/* print a sequence to a file. */
bool zSeqPrintFile(zSeq *seq, char filename[])
{
  char fname[BUFSIZ];
  FILE *fp;

  zAddSuffix( filename, (char *)ZSEQ_SUFFIX, fname, BUFSIZ );
  if( !( fp = fopen( fname, "w" ) ) ){
    ZOPENERROR( fname );
    return false;
  }
  zSeqFPrint( fp, seq );
  fclose( fp );
  return true;
}

/* print a sequence to a file. */
void zSeqFPrint(FILE *fp, zSeq *seq)
{
  zSeqCell *cp;

  zListForEachRew( seq, cp ){
    fprintf( fp, "%.10g ", cp->data.dt );
    zVecFPrint( fp, cp->data.v );
  }
}
