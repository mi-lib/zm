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
/* CLASS: zCVec
 * double precision floating point value vector class
 * ********************************************************** */

typedef struct{
  int size;
  zComplex *elem;
} zCVecStruct;
typedef zCVecStruct * zCVec;

#define _zCVecSize(v)           (v)->size
#define zCVecSize(v)            ( (v) ? _zCVecSize(v) : 0 )
#define zCVecSetSize(v,s)       ( _zCVecSize(v) = (s) )
#define zCVecSizeIsEqual(v1,v2) ( _zCVecSize(v1) == _zCVecSize(v2) )

/*! \brief abstract and set vector element.
 *
 * zCVecElem() returns the \a i th component of a vector \a v.
 * zCVecSetElem() sets the \a i th component of \a v for a
 * scalar value \a value.
 * \return
 * zCVecElem() and zCVecSetElem() return the \a i th component of \a v.
 */
#define zCVecElem(v,n)      ( &(v)->elem[n] )
#define zCVecSetElem(v,n,e) zComplexCopy( (e), zCVecElem(v,n) )

/* METHOD:
 * zCVecAlloc, zCVecFree, zCVecClear,
 * zCVecCopyNC, zCVecCopy, zCVecClone
 * - allocate, free, cleanup and copy vector.
 *
 * 'zCVecAlloc()' allocates memory for a new vector with 'size' components.
 *
 * 'zCVecFree()' destroys the given vector 'v', freeing
 * the memory allocated.
 *
 * 'zCVecClear()' clears the vector 'v', set all of the
 * components for zero.
 *
 * 'zCVecCopyNC()' copies the vector 'src' to the other
 * 'dest' without checking the size consistency between 'src'
 * and 'dest'.
 *
 * 'zCVecCopy()' copies the vector 'src' to the other 'dest'.
 *
 * 'zCVecClone()' creates a clone of the vector 'src'.
 * 's' is the size of the array.
 * \return
 * 'zCVecAlloc()' and 'zCVecClone()' return a pointer to the newly
 * created vector.
 *
 * 'zCVecFree()' returns no values.
 *
 * 'zCVecClear()' returns a pointer 'v'.
 *
 * 'zCVecCopyNC()' returns a pointer 'dest'.
 *
 * 'zCVecCopy()' returns a pointer 'dest', or the null
 * pointer if the size of 'src' and 'dest' do not coincide.
 * \notes
 * Since 'zCVecCopyNC()' does not check the size
 * consistency, when the size of 'src' and 'dest' are
 * different from each other, anything might happen.
 * If it is not urgent and you are not hasty, you should
 * use 'zCVecCopy()' for safety.
 */
__EXPORT zCVec zCVecAlloc(int s);
__EXPORT void zCVecFree(zCVec v);
__EXPORT zCVec zCVecClear(zCVec v);
__EXPORT zCVec zCVecCopyNC(zCVec src, zCVec dest);
__EXPORT zCVec zCVecCopy(zCVec src, zCVec dest);
__EXPORT zCVec zCVecClone(zCVec src);

__EXPORT zCVec zVec2CVec(zVec v, zCVec cv);

/* METHOD:
 * zCVecIsEqual - compare two vectors.
 *
 * 'zCVecIsEqual()' sees if the given two vector 'v1' and
 * 'v2' are equal to each other.
 * \return
 * 'zCVecIsEqual()' returns the true value if 'v1' equals
 * to 'v2', or the false value otherwise.
 */
__EXPORT bool zCVecIsEqual(zCVec v1, zCVec v2);

/* METHOD:
 * zCVecIsTol, zCVecIsTiny
 * - test if tiny vector.
 *
 * 'zCVecIsTol()' returns the true value if all the
 * components of the vector 'v' are less than 'tol', or the
 * false value otherwise.
 * 'zCVecIsTiny()' is the same with 'zCVecIsTol()'
 * except it compares each component with zTOL(defined in
 * 'zeda_misc.h') instead of 'tol'.
 * \return
 * 'zCVecIsTol()' and 'zCVecIsTiny()' return results as a boolean value.
 */
__EXPORT bool zCVecIsTol(zCVec v, double tol);
#define zCVecIsTiny(v) zCVecIsTol( v, zTOL )

/* METHOD:
 * zCVecAddNC, zCVecSubNC, zCVecRevNC, zCVecMulNC, zCVecDivNC,
 * zCVecAddNCDRC, zCVecSubNCDRC, zCVecRevNCDRC, zCVecMulNCDRC, zCVecDivNCDRC,
 * zCVecAdd, zCVecSub, zCVecRev, zCVecMul, zCVecDiv,
 * zCVecAddDRC, zCVecSubDRC, zCVecRevDRC, zCVecMulDRC, zCVecDivDRC,
 * - basic arithmetics for vector.
 *
 * 'zCVecAddNC()' and 'zCVecAdd()' add the two vectors,
 * 'v1' and 'v2', and put the result into 'v'.
 *
 * 'zCVecSubNC()' and 'zCVecSub()' subtract 'v2' from
 * 'v1', and put the result into 'v'.
 *
 * 'zCVecRevNC()' and 'zCVecRev()' reverse 'v1', and put
 * the result into 'v'.
 *
 * 'zCVecMulNC()' and 'zCVecMul()' multiply 'v1' by a
 * scalar value 'k', and put the result into 'v'.
 *
 * 'zCVecDivNC()' and 'zCVecDiv()' divide 'v1' by 'k',
 * and put the result into 'v'.
 *
 * 'zCVecAddNCDRC()' and 'zCVecAddDRC()' directly add
 * 'v2' to 'v1'.
 *
 * 'zCVecSubNCDRC()' and 'zCVecSubDRC()' directly
 * subtract 'v2' from 'v1'.
 *
 * 'zCVecRevNCDRC()' and 'zCVecRevDRC()' directly reverse
 * 'v'.
 *
 * 'zCVecMulNCDRC()' and 'zCVecMulDRC()' directly
 * multiply 'v' by 'k'.
 *
 * 'zCVecDivNCDRC()' and 'zCVecDivDRC()' directly divide
 * 'v' by 'k'.
 * \return
 * Each of all these functions returns a pointer to
 * the result.
 * \notes
 * The type of NC functions do not check the size
 * consistency. If it is not urgent and you are not
 * hasty, you should not use them.
 */
__EXPORT zCVec zCVecAddNC(zCVec v1, zCVec v2, zCVec v);
__EXPORT zCVec zCVecSubNC(zCVec v1, zCVec v2, zCVec v);
__EXPORT zCVec zCVecRevNC(zCVec v1, zCVec v);
__EXPORT zCVec zCVecMulNC(zCVec v1, zComplex *z, zCVec v);
__EXPORT zCVec zCVecDivNC(zCVec v1, zComplex *z, zCVec v);
__EXPORT zCVec zCVecCatNC(zCVec v1, zComplex *z, zCVec v2, zCVec v);

#define zCVecAddNCDRC(v1,v2)   zCVecAddNC( v1, v2, v1 )
#define zCVecSubNCDRC(v1,v2)   zCVecSubNC( v1, v2, v1 )
#define zCVecRevNCDRC(v)       zCVecRevNC( v, v )
#define zCVecMulNCDRC(v,z)     zCVecMulNC( v, z, v )
#define zCVecDivNCDRC(v,z)     zCVecDivNC( v, z, v )
#define zCVecCatNCDRC(v1,z,v2) zCVecCatNC( v1, z, v2, v1 )

__EXPORT zCVec zCVecAdd(zCVec v1, zCVec v2, zCVec v);
__EXPORT zCVec zCVecSub(zCVec v1, zCVec v2, zCVec v);
__EXPORT zCVec zCVecRev(zCVec v1, zCVec v);
__EXPORT zCVec zCVecMul(zCVec v1, zComplex *z, zCVec v);
__EXPORT zCVec zCVecDiv(zCVec v1, zComplex *z, zCVec v);
__EXPORT zCVec zCVecCat(zCVec v1, zComplex *z, zCVec v2, zCVec v);

#define zCVecAddDRC(v1,v2)   zCVecAdd( v1, v2, v1 )
#define zCVecSubDRC(v1,v2)   zCVecSub( v1, v2, v1 )
#define zCVecRevDRC(v)       zCVecRev( v, v )
#define zCVecMulDRC(v,z)     zCVecMul( v, z, v )
#define zCVecDivDRC(v,z)     zCVecDiv( v, z, v )
#define zCVecCatDRC(v1,z,v2) zCVecCat( v1, z, v2, v1 )

/* METHOD:
 * zCVecInnerProdNC, zCVecInnerProd
 * - inner product of vector.
 *
 * 'zCVecInnerProdNC()' and 'zCVecInnerProd()'
 * calculates the inner products of the two vector,
 * 'v1' and 'v2'.
 * \return
 * 'zCVecInnerProdNC()' and 'zCVecInnerProd()'
 * return the inner products calculated.
 * \notes
 * 'zCVecInnerProdNC()' does not check the size
 * consistency between 'v1' and 'v2'. If it is not urgent
 * and you are not hasty, you should use 'zCVecInnerProd()'.
 */
__EXPORT zComplex *zCVecInnerProdNC(zCVec v1, zCVec v2, zComplex *z);
__EXPORT zComplex *zCVecInnerProd(zCVec v1, zCVec v2, zComplex *z);

/* METHOD:
 * zCVecSqrNorm, zCVecNorm, zCVecNormalize, zCVecNormalizeDRC,
 * - normalize vector.
 *
 * 'zCVecSqrNorm()' calculates the squared norm of the
 * vector 'v'. And 'zCVecNorm()' calculates the norm
 * of 'v'.
 *
 * 'zCVecNormalize()' normalizes the vector 'src',
 * dividing by the norm of itself, and put the result
 * into 'dest'.
 * 'zCVecNormalizeDRC()' directly normalizes the
 * vector 'v'.
 * \return
 * 'zCVecSqrNorm()' and 'zCVecNorm()' return the value calculated.
 *
 * Each of 'zCVecNormalize()' and 'zCVecNormalizeDRC()'
 * returns the pointer to the result.
 */
__EXPORT double zCVecSqrNorm(zCVec v);
#define zCVecNorm(v)         sqrt( zCVecSqrNorm(v) )
__EXPORT zCVec zCVecNormalize(zCVec src, zCVec dest);
#define zCVecNormalizeDRC(v) zCVecNormalize(v,v)

/* METHOD:
 * zCVecFWrite, zCVecWrite
 * - output vector.
 *
 * 'zCVecFWrite()' writes the components of the given vector
 * 'v' to the current position of the file 'fp' in the following
 * format.
 *  n ( z1 z2 z3 ... zn )
 * 'zCVecWrite()' writes the components of 'v' in the same
 * format simply to the standard output.
 * [RETURN VALUE]
 * 'zCVecFWrite()' and 'zCVecWrite()' return no values.
 */
__EXPORT void zCVecFWrite(FILE *fp, zCVec v);
#define zCVecWrite(v) zCVecFWrite( stdout, v )

__END_DECLS

#endif /* __ZM_CVEC_H__ */
