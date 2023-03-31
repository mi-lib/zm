/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_list - vector list class.
 */

#ifndef __ZM_VEC_LIST_H__
#define __ZM_VEC_LIST_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zVecList
 * vector list class.
 * ********************************************************** */

zListClass( zVecList, zVecListCell, zVec );

/*! \brief insert a vector to a list at the head.
 *
 * zVecListInsertHead() inserts a clone of a given vector \a v to the head
 * of a list \a list.
 * \return
 * zVecListInsertHead() returns a pointer to the newly inserted list cell, if
 * it succeeds. Otherwise, the null pointer is returned.
 */
__EXPORT zVecListCell *zVecListInsertHead(zVecList *list, zVec v);

/*! \brief insert a vector to a list at the tail.
 *
 * zVecListInsertTail() inserts a clone of a given vector \a v to the tail
 * of a list \a list.
 * \return
 * zVecListInsertTail() returns a pointer to the newly inserted list cell, if
 * it succeeds. Otherwise, the null pointer is returned.
 */
__EXPORT zVecListCell *zVecListInsertTail(zVecList *list, zVec v);

/*! \brief destroy a vector list.
 *
 * zVecListDestroy() destroys a vector list \a list.
 */
__EXPORT void zVecListDestroy(zVecList *list);

/*! \brief nearest neighbor cell in a list to a vector.
 *
 * zVecListNN() finds the nearest neighbor to a vector \a v in a list of
 * vectors \a list. If \a dmin is not the null pointer, the distance
 * between \a v to the nearest neighbor is stored where is pointed by
 * \a dmin.
 * \return
 * zVecListNN() returns a pointer to the nearest neighbor vector.
 */
__EXPORT zVec zVecListNN(zVecList *list, zVec v, double *dmin);

/*! \brief scan/print a list of vectors from/to a file.
 *
 * zVecListFScan() scans vectors from a file stream \a fp, and creates
 * a list of them \a list.
 * zVecListScan() scans vectors from the standard input.
 *
 * zVecListFPrint() prints out information about a list of vectors \a list
 * to the current position of a stream \a fp.
 * zVecListPrint() prints a list \a list out to the standard output.
 * \return
 * zVecListFScan() and zVecListScan() return a pointer \a list if they
 * succeed. Otherwise, the null pointer is returned.
 *
 * zVecListFPrint() and zVecListPrint() return no value.
 */
__EXPORT zVecList *zVecListFScan(FILE *fp, zVecList *list);
#define zVecListScan(l) zVecListFScan( stdin, l )
__EXPORT void zVecListFPrint(FILE *fp, zVecList *list);
#define zVecListPrint(l) zVecListFPrint( stdout, l )

/* ********************************************************** */
/* CLASS: zVecAddrList
 * vector address list class.
 * ********************************************************** */

typedef zVecList zVecAddrList;
typedef zVecListCell zVecAddrListCell;

/*! \brief insert a pointer to a vector to a list at the head.
 */
__EXPORT zVecListCell *zVecAddrListInsertHead(zVecAddrList *list, zVec v);

/*! \brief insert a pointer to a vector to a list at the tail.
 */
__EXPORT zVecListCell *zVecAddrListInsertTail(zVecAddrList *list, zVec v);

/*! \brief destroy a vector list.
 */
__EXPORT void zVecAddrListDestroy(zVecList *list);

__END_DECLS

#endif /* __ZM_VEC_LIST_H__ */
