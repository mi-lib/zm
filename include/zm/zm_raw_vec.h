/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_raw_vec - raw vector and matrix : vector.
 */

#ifndef __ZM_RAW_VEC_H__
#define __ZM_RAW_VEC_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief cleanup and copy raw vector.
 *
 * zRawVecClear() clears a vector \a v, setting all
 * components for zeros.
 *
 * zRawVecTouchup() touches up \a v, namely, replace
 * all components which are less than zTOL for zeros.
 *
 * zRawVecCopy() copies a vector \a src to \a dest.
 *
 * \a size is the size of vector.
 * \return
 * zRawVecClear() returns a pointer to \a v.
 *
 * zRawVecTouchup() returns no value.
 *
 * zRawVecCopy() returns a pointer to \a dest.
 */
#define zRawVecClear(v,siz)      memset( v, 0, sizeof(double)*(siz) )
__EXPORT void zRawVecTouchup(double *v, int size);
#define zRawVecCopy(src,dst,siz) memcpy( dst, src, sizeof(double)*(siz) )

/*! \brief get/put a part of raw vector.
 *
 * zRawVecGet() partly copies \a src to \a dest from
 * the \a pos'th component.
 * zRawVecPut() copies \a dest to a part of \a src from
 * the \a pos'th component.
 * For both functions, the number of components copied
 * is \a size.
 * \return
 * zRawVecGet() and zRawVecPut() return a pointer to \a dest.
 */
#define zRawVecGet(src,pos,dst,siz) zRawVecCopy( (src)+(pos), dst, siz )
#define zRawVecPut(dst,pos,src,siz) zRawVecCopy( src, (dst)+(pos), siz )

/* METHOD:
 * zRawVecSetAll, zRawVecLinSpace, zRawVecRand, zRawVecShift
 * - create uniform vector, linear space vector and random vector.
 *
 * 'zRawVecSetAll()' sets all the components of 'v' for 'val'.
 * #
 * 'zRawVecLinSpace()' divides the range 'from'-'to' into
 * the values(the number of values is the same with the size
 * of 'v'), and set them as components of 'v'.
 * Each space between the components are the same.
 * #
 * 'zRawVecRand()' sets all the components of 'v' randomly
 * within the range from 'min' to 'max'.
 * #
 * 'zRawVecShift()' shifts 'v' by a scalar constant 'shift',
 * namely, 'shift' is added to all the components of 'v'.
 * #
 * 'size' is the size of vector.
 * [RETURN VALUE]
 * 'zRawVecSetAll()', 'zRawVecLinSpace()', 'zRawVecRand()'
 * and 'zRawVecShift()' return no values.
 */
__EXPORT void zRawVecSetAll(double *v, int size, double val);
__EXPORT void zRawVecLinSpace(double *v, int size, double from, double to);
__EXPORT void zRawVecRandUniform(double *v, int size, double min, double max);
__EXPORT void zRawVecRand(double *v, double *min, double *max, int size);
__EXPORT void zRawVecShift(double *v, int size, double shift);

/* METHOD:
 * zRawVecSwap - swap raw vector components.
 *
 * 'zRawVecSwap()' swaps 'i1'th value and 'i2'th value
 * of 'v'.
 * 'size' is the size of vector.
 * [RETURN VALUE]
 * 'zRawVecSwap()' returns a pointer 'v'.
 */
__EXPORT double *zRawVecSwap(double *v, int i1, int i2);

/* METHOD:
 * zRawVecIsTol, zRawVecIsTiny
 * - test for tiny raw vector.
 *
 * 'zRawVecIsTol()' returns the true value if every
 * components of 'v' are smaller than 'tol', or the
 * false value, otherwise.
 * 'zRawVecIsTiny()' compares each component of 'v'
 * with zTOL(defined in "zm_misc.h"), and returns
 * the same result with 'zRawVecIsTol()'.
 * 'size' is the size of vector.
 * [RETURN VALUE]
 * 'zRawVecIsTol()', 'zRawVecIsTiny()' return results
 * as boolean values.
 */
__EXPORT bool zRawVecIsTol(double *v, int size, double tol);
#define zRawVecIsTiny(v,siz)  zRawVecIsTol( v, siz, zTOL )

__EXPORT bool zRawVecIsNan(double *v, int size);

/* METHOD:
 * zRawVecAdd, zRawVecSub, zRawVecRev,
 * zRawVecMul, zRawVecDiv, zRawVecAmp, zRawVecDem, zRawVecCat,
 * zRawVecAddDRC, zRawVecSubDRC, zRawVecRevDRC,
 * zRawVecMulDRC, zRawVecDivDRC, zRawVecAmpDRC, zRawVecDemDRC,
 * zRawVecCatDRC, zRawVecCats, zRawVecLS
 * - basic arithmetics of raw vector.
 *
 * 'zRawVecAdd()' adds two vectors 'v1' and 'v2'.
 * 'zRawVecSub()' subtracts 'v2' from 'v1'.
 * 'zRawVecRev()' reverses the sign of 'v1'.
 * 'zRawVecMul()' multiplies 'v1' by a scalar value 'k'.
 * 'zRawVecDiv()' divides 'v1' by 'k'.
 * 'zRawVecAmp()' amplifies 'v1' by another vector 'amp'.
 * Namely, each component of 'v1' is multiplied by the
 * corresponding component of 'amp'.
 * 'zRawVecDem()' demagnifies 'v1' by another vector 'dem'.
 * Namely, each component of 'v1' is divided by the
 * corresponding component of 'dem'.
 * 'zRawVecCat()' concatenates 'v1', adding multiplied
 * 'v2' by 'k'.
 * These functions put the result into 'v'.
 * #
 * 'zRawVecAddDRC()', 'zRawVecSubDRC()', 'zRawVecRevDRC()',
 * 'zRawVecMulDRC()', 'zRawVecDivDRC()', 'zRawVecAmpDRC()',
 * 'zRawVecDemDRC()' and 'zRawVecCatDRC()' directly modify
 * the vector given as the first argument.
 * #
 * 'zRawVecCats()' concatenates 'n' vectors directly
 * to 'v'. Arguments follow 'n' as:
 *   'k1', 'v1', 'k2', 'v2', ...
 * where 'kx's are scalar values and 'vx' are vectors.
 * Then, the resultant 'v' will be:
 *   'v' + 'k1'*'v1' + 'k2'*'v2' + ...
 * 'zRawVecLS()' computes linear summation of vectors.
 * Arguments follow 'n' as:
 *   'k1', 'v1', 'k2', 'v2', ...
 * The resultant 'v' will be:
 *   'k1'*'v1' + 'k2'*'v2' + ...
 * #
 * 'size' is the size of vector.
 * [RETURN VALUE]
 * These functions return no value.
 */
__EXPORT void zRawVecAdd(double *v1, double *v2, double *v, int size);
__EXPORT void zRawVecSub(double *v1, double *v2, double *v, int size);
__EXPORT void zRawVecRev(double *v1, double *v, int size);
__EXPORT void zRawVecMul(double *v1, double k, double *v, int size);
__EXPORT void zRawVecDiv(double *v1, double k, double *v, int size);
__EXPORT void zRawVecAmp(double *v1, double *amp, double *v, int size);
__EXPORT void zRawVecDem(double *v1, double *dem, double *v, int size);
__EXPORT void zRawVecCat(double *v1, double k, double *v2, double *v, int size);
__EXPORT void zRawVecAddDRC(double *v1, double *v2, int size);
__EXPORT void zRawVecSubDRC(double *v1, double *v2, int size);
__EXPORT void zRawVecRevDRC(double *v, int size);
__EXPORT void zRawVecMulDRC(double *v, double k, int size);
__EXPORT void zRawVecDivDRC(double *v, double k, int size);
__EXPORT void zRawVecAmpDRC(double *v, double *amp, int size);
__EXPORT void zRawVecDemDRC(double *v, double *dem, int size);
__EXPORT void zRawVecCatDRC(double *v1, double k, double *v2, int size);
__EXPORT void zRawVecCats(double *v, int size, int n, ...);
__EXPORT void zRawVecLS(double *v, int size, int n, ...);

/* METHOD:
 * zRawVecInnerProd
 * - inner product of raw vector.
 *
 * 'zRawVecInnerProd()' calculates the inner product of
 * two vectors 'v1' and 'v2'.
 * 'size' is the size of vector.
 * [RETURN VALUE]
 * 'zRawVecInnerProd()' returns the inner products
 * calculated.
 */
__EXPORT double zRawVecInnerProd(double *v1, double *v2, int size);

/* METHOD:
 * zRawVecSqrNorm, zRawVecNorm,
 * zRawVecNormalize, zRawVecNormalizeDRC,
 * zRawVecSqrDist, zRawVecDist
 * - normalize raw vector.
 *
 * 'zRawVecSqrNorm()' calculates the squared norm of 'v'.
 * And, 'zRawVecNorm()' calculates the norm of 'v'.
 * #
 * 'zRawVecNormalize()' normalizes the vector 'src',
 * dividing by the norm of itself, and put the result
 * into 'dest'.
 * 'zRawVecNormalizeDRC()' directly normalizes the
 * vector 'v'.
 * #
 * 'zRawVecSqrDist()' calculates the squared distance
 * between 'v1' and 'v2'. And, 'zRawVecDist()' calculates
 * the distance between the two.
 * #
 * 'size' is the size of vector.
 * [RETURN_VALUE]
 * 'zRawVecSqrNorm()', 'zRawVecNorm()', 'zRawVecSqrDist()'
 * and 'zRawVecDist()' return the value calculated.
 * #
 * 'zRawVecNormalize()' and 'zRawVecNormalizeDRC()'
 * return the pointer to the result.
 */
__EXPORT double zRawVecSqrNorm(double *v, int size);
#define zRawVecNorm(v,siz)         sqrt( zRawVecSqrNorm(v,siz) )
__EXPORT double zRawVecWSqrNorm(double *v, double *w, int size);
#define zRawVecWNorm(v,w,siz)      sqrt( zRawVecWSqrNorm(v,w,siz) )
__EXPORT double *zRawVecNormalize(double *src, int size, double *dest);
#define zRawVecNormalizeDRC(v,siz) zRawVecNormalize(v,siz,v)
__EXPORT double zRawVecSqrDist(double *v1, double *v2, int size);
#define zRawVecDist(v1,v2,siz)     sqrt( zRawVecSqrDist( v1, v2, siz ) )

/* METHOD:
 * zRawVecFWrite, zRawVecWrite
 * - output raw vector.
 *
 * 'zRawVecFWrite()' outputs all components of a given raw
 * array of floating-point values pointed  by 'v' with a size
 * 'size' to the current position of file 'fp'.
 * 'zRawVecWrite()' outputs all components of 'v' simply to
 * the standard output.
 * [RETURN VALUE]
 * Neither 'zRawVecFWrite()' nor 'zRawVecWrite()' return
 * any values.
 */
__EXPORT void zRawVecFWrite(FILE *fp, double *v, int size);
#define zRawVecWrite(v,s) zRawVecFWrite( stdout, v, s )

__END_DECLS

#endif /* __ZM_RAW_VEC_H__ */
