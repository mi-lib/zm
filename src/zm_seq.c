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

/* zSeqInit
 * - initialize sequence.
 */
zSeq *zSeqInit(zSeq *seq)
{
  zListInit( seq );
  zListRoot(seq)->data.dt = 0;
  zListRoot(seq)->data.v = NULL;
  return seq;
}

/* zSeqFree
 * - free sequence.
 */
void zSeqFree(zSeq *seq)
{
  zSeqListCell *cp;

  while( !zListIsEmpty( seq ) )
    if( ( cp = zSeqDequeue( seq ) ) )
      zSeqListCellFree( cp );
}

/* zSeqEnqueue
 * - enqueue of a vector.
 */
zSeqListCell *zSeqEnqueue(zSeq *seq, zVec v, double dt)
{
  zSeqListCell *cp;

  if( !v ){
    ZRUNERROR( ZM_ERR_NULLVEC );
    return NULL;
  }
  if( !( cp = zAlloc( zSeqListCell, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  cp->data.v = v;
  cp->data.dt = dt;
  zQueueEnqueue( seq, cp );
  return cp;
}

/* zSeqDequeue
 * - dequeue of a vector.
 */
zSeqListCell *zSeqDequeue(zSeq *seq)
{
  zSeqListCell *cp;

  if( zListIsEmpty( seq ) ){
    ZRUNWARN( ZM_WARN_SEQ_EMPTY );
    return NULL;
  }
  zQueueDequeue( seq, &cp );
  return cp;
}

/* zSeqJump
 * - jump to the frame of a sequence.
 */
zSeqListCell *zSeqJump(zSeq *seq, int step)
{
  zSeqListCell *cp;

  for( cp=zListHead(seq); ; step--, cp=zListCellPrev(cp) ){
    if( cp == zListRoot(seq) ){
      ZRUNWARN( ZM_WARN_SEQ_STEP );
      return zListTail(seq);
    }
    if( step <= 0 ) break;
  }
  return cp;
}

/* zSeqReadFile
 * - input of sequence.
 */
bool zSeqReadFile(zSeq *seq, char filename[])
{
  FILE *fp;

  if( !( fp = zOpenFile( filename, ZSEQ_SUFFIX, "r" ) ) )
    return false;
  zSeqFRead( fp, seq );
  fclose( fp );
  return true;
}

/* zSeqFRead
 * - input of sequence.
 */
zSeq *zSeqFRead(FILE *fp, zSeq *seq)
{
  zSeqListCell *cp;

  zSeqInit( seq );
  while( !feof(fp) ){
    if( !zFSkipDelimiter( fp ) ) break;
    if( !( cp = zAlloc( zSeqListCell, 1 ) ) ){
      ZALLOCERROR();
      break;
    }
    cp->data.dt = zFDouble(fp);
    if( !( cp->data.v = zVecFRead(fp) ) ){
      ZALLOCERROR();
      free( cp );
      break;
    }
    zQueueEnqueue( seq, cp );
  }
  return seq;
}

/* zSeqWriteFile
 * - output of sequence.
 */
bool zSeqWriteFile(zSeq *seq, char filename[])
{
  char fname[BUFSIZ];
  FILE *fp;

  zAddSuffix( filename, ZSEQ_SUFFIX, fname, BUFSIZ );
  if( !( fp = fopen( fname, "w" ) ) ){
    ZOPENERROR( fname );
    return false;
  }
  zSeqFWrite( fp, seq );
  fclose( fp );
  return true;
}

/* zSeqFWrite
 * - output of sequence on .zvs format.
 */
void zSeqFWrite(FILE *fp, zSeq *seq)
{
  zSeqListCell *cp;

  zListForEachRew( seq, cp ){
    fprintf( fp, "%.10g ", cp->data.dt );
    zVecFWrite( fp, cp->data.v );
  }
}
