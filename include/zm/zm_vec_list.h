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
__ZM_EXPORT zVecListCell *zVecListInsertHead(zVecList *list, zVec v);

/*! \brief insert a vector to a list at the tail.
 *
 * zVecListInsertTail() inserts a clone of a given vector \a v to the tail
 * of a list \a list.
 * \return
 * zVecListInsertTail() returns a pointer to the newly inserted list cell, if
 * it succeeds. Otherwise, the null pointer is returned.
 */
__ZM_EXPORT zVecListCell *zVecListInsertTail(zVecList *list, zVec v);

/*! \brief destroy a vector list.
 *
 * zVecListDestroy() destroys a vector list \a list.
 */
__ZM_EXPORT void zVecListDestroy(zVecList *list);

/*! \brief select a vector of a vector list randomly.
 *
 * zVecListSelectRand() selects a vector of a vector list \a vl randomly.
 * \return
 * zVecListSelectRand() returns a pointer to the selected vector.
 */
__ZM_EXPORT zVec zVecListSelectRand(zVecList *vl);

/*! \brief nearest neighbor cell in a list to a vector.
 *
 * zVecListNN() finds the nearest neighbor to a vector \a v in a list of
 * vectors \a list. The pointer to the nearest neighbor vector is stored in \a nn.
 * \return
 * zVecListNN() returns the distance between \a v and the nearest neighbor.
 */
__ZM_EXPORT double zVecListNN(const zVecList *list, const zVec v, zVec *nn);

/* Ramer-Douglas-Peucker algorithm.
 *
 * zVecListRDP() directly thins-out elements of a list of vectors \a list based on the Ramer-Douglas-Peucker
 * algorithm, where \a metric is a pointer to a metric function that computes the distance from a sample to
 * a baseline edge, \a util is a pointer to utility data to compute the distance, and \a tol is the tolerance.
 * Namely, an element of \a list is thinned-out if the distance from it to the baseline edge that connects
 * adjacent elements is smaller than \a tol.
 * \return
 * zVecListRDP() does not return any value.
 */
__ZM_EXPORT void zVecListRDP(zVecList *list, double (* metric)(const zVec,const zVec,const zVec,void*), void *util, double tol);

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
__ZM_EXPORT zVecList *zVecListFScan(FILE *fp, zVecList *list);
#define zVecListScan(l) zVecListFScan( stdin, l )
__ZM_EXPORT void zVecListFPrint(FILE *fp, const zVecList *list);
#define zVecListPrint(l) zVecListFPrint( stdout, l )

/* ********************************************************** */
/* CLASS: zVecAddrList
 * vector address list class.
 * ********************************************************** */

typedef zVecList zVecAddrList;
typedef zVecListCell zVecAddrListCell;

/*! \brief insert a pointer to a vector to a list at the head.
 */
__ZM_EXPORT zVecListCell *zVecAddrListInsertHead(zVecAddrList *list, zVec v);

/*! \brief insert a pointer to a vector to a list at the tail.
 */
__ZM_EXPORT zVecListCell *zVecAddrListInsertTail(zVecAddrList *list, zVec v);

/*! \brief destroy a vector list.
 */
__ZM_EXPORT void zVecAddrListDestroy(zVecList *list);

__END_DECLS

#endif /* __ZM_VEC_LIST_H__ */
