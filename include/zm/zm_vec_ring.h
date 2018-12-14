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

/* METHOD:
 * zVecRingAlloc, zVecRingFree
 * - create and destroy vector ring.
 *
 * 'zVecRingAlloc()' creates a vector ring 'ring'
 * with a vector size 'dim' and a buffer size 'num'.
 * #
 * 'zVecRingFree()' destroys 'ring'.
 * [RETURN VALUE]
 * 'zVecRingAlloc()' returns the true value if succeed,
 * or the false value otherwise.
 * #
 * 'zVecRingFree()' returns no value.
 */
__EXPORT bool zVecRingAlloc(zVecRing *ring, int dim, int num);
__EXPORT void zVecRingFree(zVecRing *ring);

/* zVecRingFill
 * - fill the ring with the same vector.
 *
 * 'zVecRingFill()' fills the vector ring 'ring' monotony with
 * a vector 'v', namely, sets all the vectors of 'ring' for 'v'.
 * [RETURN VALUE]
 * 'zVecRingFill()' returns a pointer 'v' if succeeds.
 * If the size of 'v' and the components of 'ring' are different
 * from each other, the null pointer is returned.
 */
__EXPORT zVec zVecRingFill(zVecRing *ring, zVec v);

/* METHOD:
 * zVecRingCat, zVecRingLS
 * - concatenate vector ring.
 *
 * 'zVecRingCat()' concatenates vectors in the ring 'ring'
 * multiplied by correspoinding elements of 'c' directly
 * to the vector 'v', namely:
 *   'v' = 'v' + 'c1'*'ring1' + 'c2'*'ring2' + ...
 * The number of vectors summed up is determined by the size
 * of 'c'.
 * #
 * 'zVecRingLS()' computes the linear sum of 'ring', multiplied
 * by correspoinding elements of 'c'.
 * The result is put into 'v', namely:
 *   'v' = 'c1'*'ring1' + 'c2'*'ring2' + ...
 * [RETURN VALUE]
 * 'zVecRingCat()' and 'zVecRingLS()' return the resultant
 * vector 'v'.
 */
__EXPORT zVec zVecRingCat(zVec v, zVec c, zVecRing *ring);
__EXPORT zVec zVecRingLS(zVec v, zVec c, zVecRing *ring);

__END_DECLS

#endif /* __ZM_VEC_RING_H__ */
