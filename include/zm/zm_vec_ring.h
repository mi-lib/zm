/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_ring - vector ring class.
 */

#ifndef __ZM_VEC_RING_H__
#define __ZM_VEC_RING_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zVecRing
 * vector ring class.
 * ********************************************************** */

zRingClass( zVecRing, zVec );

/* create and destroy a vector ring.
 *
 * zVecRingAlloc() creates a vector ring \a ring with \a dim
 * dimensional vector. \a size is the size of the buffer.
 *
 * zVecRingFree() destroys \a ring.
 * \return
 * zVecRingAlloc() returns the true value if succeed,
 * or the false value otherwise.
 *
 * zVecRingFree() returns no value.
 */
__EXPORT bool zVecRingAlloc(zVecRing *ring, int dim, int size);
__EXPORT void zVecRingFree(zVecRing *ring);

/* fill the ring with the same vector.
 *
 * zVecRingFill() fills the vector ring \a ring monotony with
 * a vector \a v, namely, sets all the vectors of \a ring for \a v.
 * \return
 * zVecRingFill() returns a pointer \a v if succeeds.
 * If the size of \a v and the components of \a ring are different
 * from each other, the null pointer is returned.
 */
__EXPORT zVec zVecRingFill(zVecRing *ring, zVec v);

/* concatenate vector ring.
 *
 * zVecRingCat() concatenates vectors in the ring \a ring
 * multiplied by correspoinding elements of \a c directly
 * to the vector \a v, namely:
 *   \a v = \a v + \a c1 * \a ring1 + \a c2 * \a ring2 + ...
 * The number of vectors summed up is determined by the size
 * of \a c.
 *
 * zVecRingLS() computes the linear sum of \a ring, multiplied
 * by correspoinding elements of \a c.
 * The result is put into \a v, namely:
 *   \a v = \a c1 * \a ring1 + \a c2 * \a ring2 + ...
 * \return
 * zVecRingCat() and zVecRingLS() return the resultant vector
 * \a v.
 */
__EXPORT zVec zVecRingCat(zVec v, zVec c, zVecRing *ring);
__EXPORT zVec zVecRingLS(zVec v, zVec c, zVecRing *ring);

__END_DECLS

#endif /* __ZM_VEC_RING_H__ */
