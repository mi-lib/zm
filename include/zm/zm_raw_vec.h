/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_raw_vec - raw vector and matrix : vector.
 */

#ifndef __ZM_RAW_VEC_H__
#define __ZM_RAW_VEC_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief zero a raw vector.
 *
 * zRawVecZero() sets all components of a raw vector \a v for
 * zeros.
 *
 * zRawVecTouchup() replaces all components less than zTOL of
 * a raw vector \a v for zeros.
 *
 * \a size is the size of the vector.
 * \return
 * zRawVecZero() returns a pointer \a v.
 * zRawVecTouchup() returns no value.
 */
#define zRawVecZero(v,s) memset( v, 0, sizeof(double)*(s) )
__ZM_EXPORT void zRawVecTouchup(double *v, int size);

/*! \brief copy a raw vector.
 *
 * zRawVecCopy() copies a vector \a src to \a dest.
 * \a s is the size of the vectors.
 * \return
 * zRawVecCopy() returns a pointer to \a dest.
 */
#define zRawVecCopy(src,dst,s) memcpy( dst, src, sizeof(double)*(s) )

/*! \brief get/put a part of a raw vector.
 *
 * zRawVecGet() partly copies \a src to \a dest from
 * the \a pos'th component.
 * zRawVecPut() copies \a dest to a part of \a src from
 * the \a pos'th component.
 * For both functions, \a size is the number of components
 * to be copied.
 * \return
 * zRawVecGet() and zRawVecPut() return a pointer to \a dest.
 */
#define zRawVecGet(src,pos,dst,siz) zRawVecCopy( (src)+(pos), dst, siz )
#define zRawVecPut(dst,pos,src,siz) zRawVecCopy( src, (dst)+(pos), siz )

/*! \brief create uniform vector, linear space vector and random vector.
 *
 * zRawVecSetAll() sets all components of a raw vector \a v
 * for \a val.
 *
 * zRawVecLinSpace() divides the range \a from-\a to into
 * the values (the number of values is the same with the size
 * of \a v), and sets them as components of \a v.
 * Each space between the components are the same.
 *
 * zRawVecRand() sets all components of \a v randomly
 * within the range from \a min to \a max.
 *
 * \a size is the size of the raw vector.
 * \return
 * zRawVecSetAll(), zRawVecLinSpace() and zRawVecRand() return
 * no values.
 */
__ZM_EXPORT void zRawVecSetAll(double *v, int size, double val);
__ZM_EXPORT void zRawVecLinSpace(double *v, int size, double from, double to);
__ZM_EXPORT void zRawVecRandUniform(double *v, int size, double min, double max);
__ZM_EXPORT void zRawVecRand(double *v, double *min, double *max, int size);

/*! \brief shift a raw vector by a scalar value.
 *
 * zRawVecShift() shifts all elements of a raw vector \a v
 * by a constant scalar value \a shift.
 * \a size is the size of the raw vector.
 * \return
 * zRawVecShift() returns no value.
 */
__ZM_EXPORT void zRawVecShift(double *v, int size, double shift);

/*! \brief swap raw vector components.
 *
 * zRawVecSwap() swaps \a i1'th value and \a i2'th value
 * of a raw vector \a v.
 * \a size is the size of the raw vector.
 * \return
 * zRawVecSwap() returns a pointer \a v.
 */
__ZM_EXPORT double *zRawVecSwap(double *v, int i1, int i2);

/*! \brief check if a raw vector is tiny.
 *
 * zRawVecIsTol() returns the true value if all components
 * of a raw vector \a v are smaller than \a tol, or the
 * false value, otherwise.
 * zRawVecIsTiny() returns the true value if all components
 * of \a v are smaller than zTOL, which is defined in zm_misc.h,
 * or the false value, otherwise.
 *
 * \a size is the size of the raw vector.
 * \return
 * zRawVecIsTol(), zRawVecIsTiny() return results as
 * boolean values.
 */
__ZM_EXPORT bool zRawVecIsTol(double *v, int size, double tol);
#define zRawVecIsTiny(v,siz)  zRawVecIsTol( v, siz, zTOL )

/*! \brief check if a raw vector contains NaN.
 *
 * \return
 * zRawVecIsNan() returns the true value if at least one of
 * elements of a raw vector \a v is NaN, or the false value,
 * otherwise.
 * \a size is the size of the raw vector.
 */
__ZM_EXPORT bool zRawVecIsNan(double *v, int size);

/*! \brief basic arithmetics of raw vectors.
 *
 * zRawVecAdd() adds two vectors \a v1 and \a v2.
 * zRawVecSub() subtracts \a v2 from \a v1.
 * zRawVecRev() reverses the sign of \a v1.
 * zRawVecMul() multiplies \a v1 by a scalar value \a k.
 * zRawVecDiv() divides \a v1 by \a k.
 * zRawVecAmp() amplifies \a v1 by another vector \a amp.
 * Namely, each component of \a v1 is multiplied by the
 * corresponding component of \a amp.
 * zRawVecDem() demagnifies \a v1 by another vector \a dem.
 * Namely, each component of \a v1 is divided by the
 * corresponding component of \a dem.
 * zRawVecCat() concatenates \a v1 with \a v2 multiplied
 * by \a k.
 * These functions put the result into \a v.
 *
 * zRawVecAddDRC(), zRawVecSubDRC(), zRawVecRevDRC(),
 * zRawVecMulDRC(), zRawVecDivDRC(), zRawVecAmpDRC(),
 * zRawVecDemDRC() and zRawVecCatDRC() directly modify
 * the vector given as the first argument.
 *
 * zRawVecCats() concatenates \a n vectors directly
 * to \a v. Arguments follow \a n as:
 *   \a k1, \a v1, \a k2, \a v2, ...
 * where \a kx s are scalar values and \a vx s are vectors.
 * Then, the resultant \a v will be:
 *   \a v + \a k1*\a v1 + \a k2*\a v2 + ...
 * zRawVecLS() computes linear summation of vectors.
 * Arguments follow \a n as:
 *   \a k1, \a v1, \a k2, \a v2, ...
 * The resultant \a v will be:
 *   \a k1*\a v1 + \a k2*\a v2 + ...
 *
 * \a size is the size of vector.
 * \return
 * These functions return no value.
 */
__ZM_EXPORT void zRawVecAdd(double *v1, double *v2, double *v, int size);
__ZM_EXPORT void zRawVecSub(double *v1, double *v2, double *v, int size);
__ZM_EXPORT void zRawVecRev(double *v1, double *v, int size);
__ZM_EXPORT void zRawVecMul(double *v1, double k, double *v, int size);
__ZM_EXPORT void zRawVecDiv(double *v1, double k, double *v, int size);
__ZM_EXPORT void zRawVecAmp(double *v1, double *amp, double *v, int size);
__ZM_EXPORT void zRawVecDem(double *v1, double *dem, double *v, int size);
__ZM_EXPORT void zRawVecCat(double *v1, double k, double *v2, double *v, int size);
__ZM_EXPORT void zRawVecAddDRC(double *v1, double *v2, int size);
__ZM_EXPORT void zRawVecSubDRC(double *v1, double *v2, int size);
__ZM_EXPORT void zRawVecRevDRC(double *v, int size);
__ZM_EXPORT void zRawVecMulDRC(double *v, double k, int size);
__ZM_EXPORT void zRawVecDivDRC(double *v, double k, int size);
__ZM_EXPORT void zRawVecAmpDRC(double *v, double *amp, int size);
__ZM_EXPORT void zRawVecDemDRC(double *v, double *dem, int size);
__ZM_EXPORT void zRawVecCatDRC(double *v1, double k, double *v2, int size);
__ZM_EXPORT void zRawVecCats(double *v, int size, int n, ...);
__ZM_EXPORT void zRawVecLS(double *v, int size, int n, ...);

/*! \brief interior division of two raw vectors.
 *
 * zRawVecInterDiv() calculates the interior division of two
 * raw vectors \a v1 and \a v2 with a division ratio \a ratio.
 * The result is put into \a v.
 *
 * i.e. \a v = (1-\a ratio)* \a v1 + \a ratio * \a v2.
 * \return
 * zRawVecInterDiv() does not return any value.
 */
__ZM_EXPORT void zRawVecInterDiv(double *v1, double *v2, double ratio, double *v, int size);

/*! \brief midpoint of two raw vectors.
 *
 * zRawVecMid() calculates the midpoint of two raw vectors \a v1 and
 * \a v2. The result is put into \a v.
 *
 * i.e. \a v = ( \a v1 + \a v2 ) / 2.
 * \return
 * zRawVecMid() does not return any value.
 */
__ZM_EXPORT void zRawVecMid(double *v1, double *v2, double *v, int size);

/*! \brief scale a raw vector with two boundary vectors.
 *
 * zRawVecScale() scales a raw vector \a x with the minimum and maximum
 * boundary vectors \a min and \a max, respectively and puts it into
 * \a v. Namely, each component of \a v is the corresponding component
 * of \a x multipied by the corresponding components of \a max - \a min
 * and offset by that of \a min.
 * \a size is the size of the vectors.
 * \return
 * zRawVecScale() does not return any value.
 */
__ZM_EXPORT void zRawVecScale(double *x, double *min, double *max, double *v, int size);

/*! \brief inner product of raw vector.
 *
 * \return
 * zRawVecInnerProd() returns the inner product of
 * two raw vectors \a v1 and \a v2.
 * \a size is the size of the vectors.
 */
__ZM_EXPORT double zRawVecInnerProd(double *v1, double *v2, int size);

/*! \brief normalize a raw vector.
 *
 * zRawVecSqrNorm() returns the squared norm of a raw
 * vector \a v.
 * zRawVecNorm() returns the norm of \a v.
 *
 * zRawVecNormalize() normalizes a raw vector \a src,
 * namely, \a src is divided by the norm of itself.
 * The result is put into \a dest.
 * zRawVecNormalizeDRC() directly normalizes \a v.
 *
 * zRawVecSqrDist() returns the squared distance between
 * \a v1 and \a v2.
 * zRawVecDist() returns the distance between \a v1 and
 * \a v2.
 *
 * \a size is the size of the vector(s).
 * \return
 * zRawVecSqrNorm(), zRawVecNorm(), zRawVecSqrDist()
 * and zRawVecDist() return the value calculated.
 *
 * zRawVecNormalize() and zRawVecNormalizeDRC() return
 * a pointer to the result.
 */
__ZM_EXPORT double zRawVecSqrNorm(double *v, int size);
#define zRawVecNorm(v,siz)         sqrt( zRawVecSqrNorm(v,siz) )
__ZM_EXPORT double zRawVecWSqrNorm(double *v, double *w, int size);
#define zRawVecWNorm(v,w,siz)      sqrt( zRawVecWSqrNorm(v,w,siz) )
__ZM_EXPORT double *zRawVecNormalize(double *src, int size, double *dest);
#define zRawVecNormalizeDRC(v,siz) zRawVecNormalize(v,siz,v)
__ZM_EXPORT double zRawVecSqrDist(double *v1, double *v2, int size);
#define zRawVecDist(v1,v2,siz)     sqrt( zRawVecSqrDist( v1, v2, siz ) )

/*! \brief print a raw vector.
 *
 * zRawVecFPrint() prints all components of a raw vector
 * \a v to the current position of a file \a fp.
 * \a size is the size of the vector.
 * zRawVecPrint() prints all components of \a v out to
 * the standard output.
 * \return
 * zRawVecFPrint() and zRawVecPrint() return no values.
 */
__ZM_EXPORT void zRawVecFPrint(FILE *fp, double *v, int size);
#define zRawVecPrint(v,s) zRawVecFPrint( stdout, v, s )

__END_DECLS

#endif /* __ZM_RAW_VEC_H__ */
