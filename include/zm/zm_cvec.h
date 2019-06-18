/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_cvec - complex vector class.
 */

#ifndef __ZM_CVEC_H__
#define __ZM_CVEC_H__

#include <zm/zm_vec.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \brief complex number vector class
 * ********************************************************** */
zArrayClass( zCVecStruct, zComplex );
typedef zCVecStruct * zCVec;

/*! \brief size of a complex vector.
 * \retval the size of a complex vector if \a v is not null.
 * \retval 0 if \a v is the null pointer.
 */
#define zCVecSizeNC(v)          zArraySize(v)
#define zCVecSize(v)            ( (v) ? zCVecSizeNC(v) : 0 )
/*! \brief set the size of a complex vector. */
#define zCVecSetSizeNC(v,s)     ( zCVecSizeNC(v) = (s) )
#define zCVecSetSize(v,s)       ( (v) ? zCVecSetSizeNC(v,s) : 0 )
/*! \brief check if the sizes of two complex vectors are equal. */
#define zCVecSizeIsEqual(v1,v2) ( zCVecSizeNC(v1) == zCVecSizeNC(v2) )

/*! \brief abstract and set an element of a complex vector.
 *
 * zCVecElem() returns the \a i th component of a vector \a v.
 * zCVecSetElem() sets the \a i th component of \a v for a
 * scalar value \a value.
 * \return
 * zCVecElem() and zCVecSetElem() return the \a i th component of \a v.
 */
#define zCVecElemNC(v,n)      zArrayElemNC(v,n)
#define zCVecElem(v,n)        zArrayElem(v,n)
#define zCVecSetElemNC(v,n,e) zArraySetElemNC(v,n,e)
#define zCVecSetElem(v,n,e)   zArraySetElem(v,n,e)

/*! \brief allocate, free, zero and copy a complex vector.
 *
 * zCVecAlloc() allocates memory for a new vector with \a size.
 *
 * zCVecFree() frees memory allocated for a complex vector \a v.
 *
 * zCVecZero() sets all components of a complex vector \a v
 * for zeros.
 *
 * zCVecCopyNC() copies a complex vector \a src to another
 * \a dest without checking the size consistency between them.
 *
 * zCVecCopy() copies a complex vector \a src to another \a dest.
 *
 * zCVecClone() creates a clone of a complex vector \a src.
 * \return
 * zCVecAlloc() and zCVecClone() return a pointer to the newly
 * allocated vector.
 *
 * zCVecFree() returns no values.
 *
 * zCVecZero() returns a pointer \a v.
 *
 * zCVecCopyNC() returns a pointer \a dest.
 *
 * zCVecCopy() returns a pointer \a dest, or the null pointer
 * if the sizes of \a src and \a dest do not coincide.
 * \notes
 * Since zCVecCopyNC() does not check the size consistency,
 * anything might happen if the sizes of \a src and \a dest
 * are different. If it is not urgent and you are not hasty,
 * you should use zCVecCopy() for safety.
 */
__EXPORT zCVec zCVecAlloc(int size);
__EXPORT void zCVecFree(zCVec v);
__EXPORT zCVec zCVecZero(zCVec v);
__EXPORT zCVec zCVecCopyNC(zCVec src, zCVec dest);
__EXPORT zCVec zCVecCopy(zCVec src, zCVec dest);
__EXPORT zCVec zCVecClone(zCVec src);

__EXPORT zCVec zVec2CVec(zVec v, zCVec cv);

/*! \brief compare two complex vectors.
 *
 * zCVecIsEqual() checks if two complex vectors \a v1 and
 * \a v2 are equal.
 * \return
 * zCVecIsEqual() returns the true value if \a v1 equals
 * to \a v2, or the false value otherwise.
 */
__EXPORT bool zCVecIsEqual(zCVec v1, zCVec v2);

/*! \brief check if a complex vector is tiny.
 *
 * zCVecIsTol() returns the true value if all components
 * of a complex vector \a v are less than \a tol, or the
 * false value otherwise.
 *
 * zCVecIsTiny() is the same with zCVecIsTol() except that
 * it uses zTOL (defined in "zeda_misc.h") for the tolerance.
 * \return
 * zCVecIsTol() and zCVecIsTiny() return a boolean value.
 */
__EXPORT bool zCVecIsTol(zCVec v, double tol);
#define zCVecIsTiny(v) zCVecIsTol( v, zTOL )

/*! \brief basic arithmetics for the complex vector.
 *
 * zCVecAddNC() and zCVecAdd() add two complex vectors \a v1
 * and \a v2, and puts the result into \a v.
 *
 * zCVecSubNC() and zCVecSub() subtract \a v2 from \a v1,
 * and puts the result into \a v.
 *
 * zCVecRevNC() and zCVecRev() reverse \a v1, and puts the
 * result into \a v.
 *
 * zCVecMulNC() and zCVecMul() multiply \a v1 by a scalar
 * value \a k, and puts the result into \a v.
 *
 * zCVecDivNC() and zCVecDiv() divide \a v1 by \a k, and
 * puts the result into \a v.
 *
 * zCVecAddNCDRC() and zCVecAddDRC() directly add \a v2 to
 * \a v1.
 *
 * zCVecSubNCDRC() and zCVecSubDRC() directly subtract \a v2
 * from \a v1.
 *
 * zCVecRevNCDRC() and zCVecRevDRC() directly reverse \a v.
 *
 * zCVecMulNCDRC() and zCVecMulDRC() directly multiply \a v
 * by \a k.
 *
 * zCVecDivNCDRC() and zCVecDivDRC() directly divide \a v
 * by \a k.
 * \return
 * These functions return a pointer to the result.
 * \notes
 * The NC-type functions do not check the size consistency.
 * If it is not urgent and you are not hasty, you should not
 * use them.
 */
__EXPORT zCVec zCVecAddNC(zCVec v1, zCVec v2, zCVec v);
__EXPORT zCVec zCVecSubNC(zCVec v1, zCVec v2, zCVec v);
__EXPORT zCVec zCVecRevNC(zCVec v1, zCVec v);
__EXPORT zCVec zCVecMulNC(zCVec v1, zComplex *z, zCVec v);
__EXPORT zCVec zCVecDivNC(zCVec v1, zComplex *z, zCVec v);
__EXPORT zCVec zCVecCatNC(zCVec v1, zComplex *z, zCVec v2, zCVec v);

__EXPORT zCVec zCVecAdd(zCVec v1, zCVec v2, zCVec v);
__EXPORT zCVec zCVecSub(zCVec v1, zCVec v2, zCVec v);
__EXPORT zCVec zCVecRev(zCVec v1, zCVec v);
__EXPORT zCVec zCVecMul(zCVec v1, zComplex *z, zCVec v);
__EXPORT zCVec zCVecDiv(zCVec v1, zComplex *z, zCVec v);
__EXPORT zCVec zCVecCat(zCVec v1, zComplex *z, zCVec v2, zCVec v);

#define zCVecAddNCDRC(v1,v2)   zCVecAddNC( v1, v2, v1 )
#define zCVecSubNCDRC(v1,v2)   zCVecSubNC( v1, v2, v1 )
#define zCVecRevNCDRC(v)       zCVecRevNC( v, v )
#define zCVecMulNCDRC(v,z)     zCVecMulNC( v, z, v )
#define zCVecDivNCDRC(v,z)     zCVecDivNC( v, z, v )
#define zCVecCatNCDRC(v1,z,v2) zCVecCatNC( v1, z, v2, v1 )

#define zCVecAddDRC(v1,v2)   zCVecAdd( v1, v2, v1 )
#define zCVecSubDRC(v1,v2)   zCVecSub( v1, v2, v1 )
#define zCVecRevDRC(v)       zCVecRev( v, v )
#define zCVecMulDRC(v,z)     zCVecMul( v, z, v )
#define zCVecDivDRC(v,z)     zCVecDiv( v, z, v )
#define zCVecCatDRC(v1,z,v2) zCVecCat( v1, z, v2, v1 )

/*! \brief inner product of two complex vectors.
 *
 * \return
 * zCVecInnerProdNC() and zCVecInnerProd() calculates
 * the inner products of two complex vectors \a v1 and \a v2.
 * \notes
 * zCVecInnerProdNC() does not check the size consistency
 * between \a v1 and \a v2. If it is not urgent and you are
 * not hasty, you should use zCVecInnerProd().
 */
__EXPORT zComplex *zCVecInnerProdNC(zCVec v1, zCVec v2, zComplex *z);
__EXPORT zComplex *zCVecInnerProd(zCVec v1, zCVec v2, zComplex *z);

/*! \brief normalize vector.
 *
 * zCVecSqrNorm() calculates the squared norm of a complex
 * vector \a v. zCVecNorm() calculates the norm of \a v.
 *
 * zCVecNormalize() normalizes a complex vector \a src,
 * namely, divides it by the norm of itself, and puts the
 * result into \a dest.
 *
 * zCVecNormalizeDRC() directly normalizes a complex vector
 * \a v.
 * \return
 * zCVecSqrNorm() and zCVecNorm() return the calculated value.
 *
 * zCVecNormalize() and zCVecNormalizeDRC() return a pointer
 * to the result.
 */
__EXPORT double zCVecSqrNorm(zCVec v);
#define zCVecNorm(v)         sqrt( zCVecSqrNorm(v) )
__EXPORT zCVec zCVecNormalize(zCVec src, zCVec dest);
#define zCVecNormalizeDRC(v) zCVecNormalize(v,v)

/*! \brief print a complex vector.
 *
 * zCVecFPrint() prints components of a complex vector \a v
 * to the current position of a file \a fp in the following
 * format.
 *  n ( z1 z2 z3 ... zn )
 *
 * zCVecPrint() prints components of \a v in the same format
 * to the standard output.
 * \return
 * zCVecFPrint() and zCVecPrint() return no values.
 */
__EXPORT void zCVecFPrint(FILE *fp, zCVec v);
#define zCVecPrint(v) zCVecFPrint( stdout, v )

__END_DECLS

#endif /* __ZM_CVEC_H__ */
