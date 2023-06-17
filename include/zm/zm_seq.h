/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_seq - vector sequence.
 */

#ifndef __ZM_SEQ_H__
#define __ZM_SEQ_H__

#include <zeda/zeda.h>
#include <zm/zm_vec.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zSeq
 * motion sequence class for multiple dimension system.
 * ********************************************************** */

typedef struct{
  double dt;
  zVec v;
} zSeqCell;

zListClass( zSeq, zSeqListCell, zSeqCell );

#define zSeqListCellFree(c) do{\
  if( (c) ){\
    zVecFree( (c)->data.v );\
    free( (c) );\
    (c) = NULL;\
  }\
} while(0)

/*! \brief initialize and free a sequence.
 *
 * zSeqInit() initializes a sequence pointed by \a seq,
 * setting the dimension for \a dim and the delta time step
 * for \a dt.
 *
 * zSeqFree() frees the internal list of a sequence pointed
 * by \a seq.
 * \return
 * zSeqInit() returns the pointer \a seq.
 * zSeqFree() returns no value.
 */
__ZM_EXPORT zSeq *zSeqInit(zSeq *seq);
__ZM_EXPORT void zSeqFree(zSeq *seq);

/*! \brief add and delete a vector cell to a sequence.
 *
 * zSeqEnqueue() enqueues a new vector at the tail of \a seq,
 * while zSeqDequeue() dequeues the vector at the head of \a seq.
 * These two functions do not modify the cursor.
 * \return
 * zSeqEnqueue() and zSeqDequeue() return a pointer to \a v,
 * if succeeding, or the null vector, otherwise.
 */
__ZM_EXPORT zSeqListCell *zSeqEnqueue(zSeq *seq, zVec v, double dt);
__ZM_EXPORT zSeqListCell *zSeqDequeue(zSeq *seq);

/*! \brief forward/backward frame step and rewind a sequence.
 *
 * zSeqJump() places the cursor at \a step'th frame from the
 * beginning of \a seq.
 * \return
 * All these functions return a pointer to the current step.
 */
__ZM_EXPORT zSeqListCell *zSeqJump(zSeq *seq, int step);

/*! \brief scan and print a sequence to a file.
 *
 * These functions scan/print a sequence from/to a file.
 * The format written in ASCII charactors is as follows.
 *   DT DIM v0_1 v0_2 v0_3 ... v0_DIM
 *   ...
 *   DT DIM vN_1 vN_2 vN_3 ... vN_DIM
 *
 * zSeqScanFile() scans a file named \a filename or \a filename.zvs
 * and creates a sequence \a seq.
 * zSeqFScan() scans a sequence from the current position
 * of a file \a fp and puts it into a sequence pointed by
 * \a seq. zSeqScan() scans a sequence from the standard
 * input.
 *
 * zSeqPrintFile() prints a sequence \a seq out to a file
 * named \a filename or \a filename.zvs.
 * zSeqFPrint() prints a sequence \a seq out to the current
 * position of a file \a fp.
 * zSeqPrint() prints a sequence out to the standard output.
 * \return
 * zSeqScanFile() and zSeqPrintFile() return the true value
 * if succeeds, or the false value otherwise.
 *
 * zSeqFScan() and zSeqScan() return a pointer to \a seq if
 * succeeds, or the null pointer otherwise.
 *
 * zSeqPrintFile(), zSeqFPrint() and zSeqPrint() return no
 * values.
 */
#define ZSEQ_SUFFIX "zvs"
__ZM_EXPORT bool zSeqScanFile(zSeq *seq, char filename[]);
__ZM_EXPORT zSeq *zSeqFScan(FILE *fp, zSeq *seq);
#define zSeqScan(s)  zSeqFScan( stdin, (s) )

__ZM_EXPORT bool zSeqPrintFile(zSeq *seq, char filename[]);
__ZM_EXPORT void zSeqFPrint(FILE *fp, zSeq *seq);
#define zSeqPrint(s) zSeqFPrint( stdout, (s) )

__END_DECLS

#endif /* __ZM_SEQ_H__ */
