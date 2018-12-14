/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_array - vector array class.
 */

#ifndef __ZM_VEC_ARRAY_H__
#define __ZM_VEC_ARRAY_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zVecArray
 * vector array class.
 * ********************************************************** */

zArrayClass( zVecArray, zVec );

/*! \brief allocate vector array.
 *
 * zVecArrayAlloc() allocate \a num \a dim -dimensional vectors in \a array.
 * \retval true if succeeding to allocate.
 * \retval false if failing to allocate.
 */
__EXPORT bool zVecArrayAlloc(zVecArray *array, int dim, int num);

/*! \brief free vector array.
 *
 * zVecArrayFree() frees a vector array \a array.
 */
__EXPORT void zVecArrayFree(zVecArray *array);

/*! \brief an element of a vector in an array */
#define zVecArrayElem(a,i,j) zVecElem( *zArrayElem( a, i ), j )

/*! \brief fill a vector array uniformly with the same vector.
 *
 * zVecArrayFill() fills a vector array \a array uniformly with
 * a vector \a v i.e. it sets all the vectors of \a array for \a v.
 * \retval \a v if succeeding.
 * \retval the null pointer if the vector size of \a array is
 * inconsistent with \a v.
 */
__EXPORT zVec zVecArrayFill(zVecArray *array, zVec v);

__END_DECLS

#endif /* __ZM_VEC_ARRAY_H__ */
