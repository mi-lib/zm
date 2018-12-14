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
 * if \a flag is true, a clone of \a v is registered.
 */
__EXPORT zVecListCell *zVecListInsertHead(zVecList *list, zVec v, bool flag);

/*! \brief insert a vector to a list at the tail.
 * if \a flag is true, a clone of \a v is registered.
 */
__EXPORT zVecListCell *zVecListInsertTail(zVecList *list, zVec v, bool flag);

/*! \brief destroy a vector list.
 * if \a flag is true, a vector pointed by each cell is freed.
 */
__EXPORT void zVecListDestroy(zVecList *list, bool flag);

/*! \brief nearest neighbor cell in a list to a vector.
 * the distance from \a v to the nearest neighbor is stored where
 * is pointed by \a dmin, unless \a dmin is the null pointer.
 */
__EXPORT zVec zVecListNN(zVecList *list, zVec v, double *dmin);

/*! \brief output a vector list to a file.
 */
__EXPORT void zVecListFWrite(FILE *fp, zVecList *list);
#define zVecListWrite(l) zVecListFWrite( stdout, l )

__END_DECLS

#endif /* __ZM_VEC_LIST_H__ */
