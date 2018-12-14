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
__EXPORT zSeq *zSeqInit(zSeq *seq);
__EXPORT void zSeqFree(zSeq *seq);

/*! \brief add and delete a vector cell to a sequence.
 *
 * zSeqEnqueue() enqueues a new vector at the tail of \a seq,
 * while zSeqDequeue() dequeues the vector at the head of \a seq.
 * These two functions do not modify the cursor.
 * \return
 * zSeqEnqueue() and zSeqDequeue() return a pointer to \a v,
 * if succeeding, or the null vector, otherwise.
 */
__EXPORT zSeqListCell *zSeqEnqueue(zSeq *seq, zVec v, double dt);
__EXPORT zSeqListCell *zSeqDequeue(zSeq *seq);

/*! \brief forward/backward frame step and rewind a sequence.
 *
 * zSeqJump() places the cursor at \a step'th frame from the
 * beginning of \a seq.
 * \return
 * All these functions return a pointer to the current step.
 */
__EXPORT zSeqListCell *zSeqJump(zSeq *seq, int step);

/*! \brief input/output a sequence to a file.
 *
 * zSeqReadFile() reads the file named \a filename or \a filename.zvs
 * and creates a sequence \a seq.
 * The format written in ASCII charactors is as follows.
 *   DT DIM v0_1 v0_2 v0_3 ... v0_DIM
 *   ...
 *   DT DIM vN_1 vN_2 vN_3 ... vN_DIM
 *
 * zSeqFRead() reads the sequence from the current position
 * of the file \a fp and puts it into the sequence pointed by
 * \a seq. zSeqRead() simply reads the sequence from the
 * standard input.
 *
 * zSeqWriteFile() writes a sequence \a seq to the file named
 * \a filename or \a filename.zvs.
 * The format is the same with one for zSeqReadFile().
 *
 * zSeqFWrite() writes a sequence \a seq to the current
 * position of the file \a fp in the same rule of the format
 * with zSeqFRead().
 *
 * zSeqWrite() writes the sequence to the standard output.
 * \return
 * zSeqReadFile() and zSeqWriteFile() return the true value
 * if succeeds, or the false value otherwise.
 *
 * Other zSeqRead() family functions return a pointer to \a seq
 * if succeeds, or the null pointer otherwise. Other functions
 * of zSeqWrite() family return no values.
 */
#define ZSEQ_SUFFIX "zvs"
__EXPORT bool zSeqReadFile(zSeq *seq, char filename[]);
__EXPORT zSeq *zSeqFRead(FILE *fp, zSeq *seq);
#define zSeqRead(s)  zSeqFRead( stdin, (s) )

__EXPORT bool zSeqWriteFile(zSeq *seq, char filename[]);
__EXPORT void zSeqFWrite(FILE *fp, zSeq *seq);
#define zSeqWrite(s) zSeqFWrite( stdout, (s) )

__END_DECLS

#endif /* __ZM_SEQ_H__ */
