/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_cmat - complex matrix class.
 */

#ifndef __ZM_CMAT_H__
#define __ZM_CMAT_H__

#include <zm/zm_mat.h>
#include <zm/zm_cvec.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zCMat
 * complex matrix class
 * NOTES: each elements of matrix(size=r*c) is at (0 - r-1,0 - c-1).
 * ********************************************************** */
zArray2Class( zCMatStruct, zComplex );
typedef zCMatStruct * zCMat;

#define zCMatRowSizeNC(m)        zArray2RowSize(m)
#define zCMatRowSize(m)          ( (m) ? zCMatRowSizeNC(m) : 0 )
#define zCMatColSizeNC(m)        zArray2ColSize(m)
#define zCMatColSize(m)          ( (m) ? zCMatColSizeNC(m) : 0 )
#define zCMatSetRowSizeNC(m,r)   ( zCMatRowSizeNC(m) = (r) )
#define zCMatSetRowSize(m,r)     ( (m) ? zCMatSetRowSizeNC(m,r) : 0 )
#define zCMatSetColSizeNC(m,c)   ( zCMatColSizeNC(m) = (c) )
#define zCMatSetColSize(m,c)     ( (m) ? zCMatSetColSizeNC(m,c) : 0 )

#define zCMatSetSizeNC(m,r,c) do{\
  zCMatSetRowSizeNC(m,r);\
  zCMatSetColSizeNC(m,c);\
} while(0)
#define zCMatSetSize(m,r,c) do{\
  zCMatSetRowSize(m,r);\
  zCMatSetColSize(m,c);\
} while(0)

#define zCMatRowSizeEqual(m1,m2)    ( zCMatRowSizeNC(m1) == zCMatRowSizeNC(m2) )
#define zCMatColSizeEqual(m1,m2)    ( zCMatColSizeNC(m1) == zCMatColSizeNC(m2) )
#define zCMatSizeEqual(m1,m2)       ( zCMatRowSizeEqual(m1,m2) && zCMatColSizeEqual(m1,m2) )
#define zCMatColCVecSizeEqual(m,v)  ( zCMatColSizeNC(m) == zCVecSizeNC(v) )
#define zCMatRowCVecSizeEqual(m,v)  ( zCMatRowSizeNC(m) == zCVecSizeNC(v) )
#define zCMatRowColSizeEqual(m1,m2) ( zCMatRowSizeNC(m1) == zCMatColSizeNC(m2) )
#define zCMatColRowSizeEqual(m1,m2) zCMatRowColSizeEqual(m2,m1)
#define zCMatIsSqr(m)               zCMatRowColSizeEqual(m,m)

/*! \brief pointer to the array buffer of a complex matrix. */
#define zCMatBufNC(m)      zArrayBuf(m)
#define zCMatBuf(m)        ( (m) ? zCMatBufNC(m) : NULL )
/*! \brief pointer to the \a r th row array buffer of a complex matrix. */
#define zCMatRowBufNC(m,r) ( zCMatBufNC(m) + (r)*zCMatColSizeNC(m) )
#define zCMatRowBuf(m,r)   ( (m) ? zCMatRowBufNC(m,r) : NULL )

/*! \brief retrieve and set an element of a complex matrix.
 *
 * zCMatElem() returns the component at \a r th row and \a c th column of a complex matrix \a m.
 *
 * zCMatSetElem() sets the component at \a r th row and \a c th column of a complex matrix \a m
 * for a scalar value \a value.
 * \return
 * zCMatElem() and zCMatSetElem() return the address of the component at \a r th row and \a c th
 * column of \a m.
 */
#define zCMatElemNC(m,r,c)      zArray2ElemNC(m,r,c)
#define zCMatElem(m,r,c)        zArray2Elem(m,r,c)
#define zCMatSetElemNC(m,r,c,e) zArray2SetElemNC(m,r,c,e)
#define zCMatSetElem(m,r,c,e)   zArray2SetElem(m,r,c,e)

/*! \brief allocate, free and zero a complex matrix.
 *
 * zCMatAlloc() allocates memory for a new matrix with the size \a row times \a col.
 *
 * zCMatAllocSqr() allocates memory for a new square complex matrix with the size \a size times \a size.
 *
 * zCMatFree() frees the memory allocated for a complex matrix \a m.
 *
 * zCMatZero() sets all components of a complex matrix \a m for zeros.
 * \return
 * zCMatAlloc() and zCMatAllocSqr() return a pointer to the newly allocated memory area.
 *
 * zCMatFree() returns no value.
 *
 * zCMatZero() returns a pointer \a m.
 */
__ZM_EXPORT zCMat zCMatAlloc(int row, int col);
#define zCMatAllocSqr(s) zCMatAlloc( (s), (s) )
__ZM_EXPORT void zCMatFree(zCMat m);
__ZM_EXPORT zCMat zCMatZero(zCMat m);

/*! \brief touchup a complex matrix.
 *
 * zCVecTouchup() replaces real part or imaginary part of all components of a complex matrix \a m
 * with zero if either value relative to the other part is smaller than \a tol.
 * \return
 * zCMatTouchup() returns the pointer \a m.
 */
__ZM_EXPORT zCMat zCMatTouchup(zCMat m, double tol);

/*! \brief create a random complex matrix with a uniform range. */
__ZM_EXPORT zCMat zCMatRandUniform(zCMat m, double rmin, double rmax, double imin, double imax);
/*! \brief create a complex random matrix with range matrices. */
__ZM_EXPORT zCMat zCMatRand(zCMat m, zCMat min, zCMat max);

/*! \brief copy a complex  matrix.
 *
 * zCMatCopyNC() copies a complex matrix \a src to the other \a dest without checking the size
 * consistency between \a src and \a dest.
 *
 * zCMatCopy() copies a complex matrix \a src to the other \a dest.
 *
 * zCMatClone() creates a clone of a complex matrix \a m, namely, allocates memory and copies \a m to it.
 * \return
 * zCMatCopyNC() returns a pointer \a dest.
 *
 * zCMatCopy() returns a pointer \a dest if succeedings, or the null pointer if the size of \a src and
 * \a dest are different.
 *
 * zCMatClone() returns a pointer to the newly created matrix.
 * \notes
 * Since zCMatCopyNC() does not check the size consistency, anything might happen if the size of \a src
 * and \a dest are different. If it is not urgent and you are not hasty, you should use zCMatCopy() for
 * safety.
 */
__ZM_EXPORT zCMat zCMatCopyNC(const zCMat src, zCMat dest);
__ZM_EXPORT zCMat zCMatCopy(const zCMat src, zCMat dest);
__ZM_EXPORT zCMat zCMatClone(const zCMat src);

/*! \brief convert a real matrix to a complex matrix. */
__ZM_EXPORT zCMat zMatToCMat(const zMat m, zCMat cm);

/*! \brief abstract and put row/column vector of a complex matrix.
 *
 * zCMatRowNC() abstracts the \a row th row vector of a complex matrix \a m and puts it into a complex
 * vector \a v without checking the size.
 * zCMatRow() abstracts the \a row th row vector of \a m and puts it into \a v.
 *
 * zCMatColNC() abstracts the \a col th column vector of \a m and puts it into \a v without checking
 * the size.
 * zCMatCol() abstracts the \a col th column vector of \a m and puts it into \a v.
 *
 * zCMatPutRowNC() puts the \a row th row vector of \a m for \a v without checking the size.
 * zCMatPutRow() puts the \a row th row vector of \a m for \a v.
 *
 * zCMatPutColNC() puts the \a col th column vector of \a m for \a v without checking the size.
 * zCMatPutCol() puts the \a col th column vector of \a m for \a v.
 * \return
 * zCMatRowNC(), zCMatColNC(), zCMatRow(), and zCMatCol() return a pointer to the abstracted vector.
 *
 * zCMatPutRowNC(), zCMatPutColNC(), zCMatPutRow(), and zCMatPutCol() return a pointer \a m.
 * \notes
 * If it is not urgent and you are not hasty, you'd better not use NC functions for safety.
 */
__ZM_EXPORT zCVec zCMatGetRowNC(const zCMat m, int row, zCVec v);
__ZM_EXPORT zCVec zCMatGetColNC(const zCMat m, int col, zCVec v);
__ZM_EXPORT zCVec zCMatGetRow(const zCMat m, int row, zCVec v);
__ZM_EXPORT zCVec zCMatGetCol(const zCMat m, int col, zCVec v);
__ZM_EXPORT zCMat zCMatPutRowNC(zCMat m, int row, const zCVec v);
__ZM_EXPORT zCMat zCMatPutColNC(zCMat m, int col, const zCVec v);
__ZM_EXPORT zCMat zCMatPutRow(zCMat m, int row, const zCVec v);
__ZM_EXPORT zCMat zCMatPutCol(zCMat m, int col, const zCVec v);

/*! \brief check if a complex matrix is tiny.
 *
 * zCMatIsTol() returns the true value if all components of a complex matrix \a m are less than \a tol.
 * Otherwise, the false value is returned.
 *
 * zCMatIsTiny() is the same with zCMatIsTol() except that it compares each component with zTOL, which
 * is defined in zm_misc.h as the tolerance.
 * \return
 * zCMatIsTol() and zCMatIsTiny() return a boolean value.
 */
__ZM_EXPORT bool zCMatIsTol(const zCMat m, double tol);
#define zCMatIsTiny(m) zCMatIsTol( (m), zTOL )

/*! \brief basic arithmetics for the complex matrix.
 *
 * zCMatAddNC() and zCMatAdd() add two complex matrices \a m1 and \a m2, and put the result into \a m.
 *
 * zCMatSubNC() and zCMatSub() subtract \a m2 from \a m1, and put the result into \a m.
 *
 * zCMatRevNC() and zCMatRev() reverse \a m1, and put the result into \a m.
 *
 * zCMatMulNC() and zCMatMul() multiply \a m1 by a real scalar value \a k, and put the result into \a m.
 *
 * zCMatDivNC() and zCMatDiv() divide \a m1 by \a k, and put the result into \a m.
 *
 * zCMatCMulNC() and zCMatCMul() multiply \a m1 by a complex scalar value \a z, and put the result into \a m.
 *
 * zCMatCDivNC() and zCMatCDiv() divide \a m1 by \a z, and put the result into \a m.
 *
 * zCMatAddNCDRC() and zCMatAddDRC() directly add \a m2 to \a m1.
 *
 * zCMatSubNCDRC() and zCMatSubDRC() directly subtract \a m2 from \a m1.
 *
 * zCMatRevNCDRC() and zCMatRevDRC() directly reverse \a m.
 *
 * zCMatMulNCDRC() and zCMatMulDRC() directly multiply \a m by \a k.
 *
 * zCMatDivNCDRC() and zCMatDivDRC() directly divide \a m by \a k.
 *
 * zCMatCMulNCDRC() and zCMatCMulDRC() directly multiply \a m by \a z.
 *
 * zCMatCDivNCDRC() and zCMatCDivDRC() directly divide \a m by \a z.
 * \return
 * These functions return a pointer to the result.
 * \notes
 * NC-type functions do calculation without checking the size consistency. If it is not urgent and
 * you are not hasty, you should not use them.
 */
__ZM_EXPORT zCMat zCMatAddNC(const zCMat m1, const zCMat m2, zCMat m);
__ZM_EXPORT zCMat zCMatSubNC(const zCMat m1, const zCMat m2, zCMat m);
__ZM_EXPORT zCMat zCMatRevNC(const zCMat m1, zCMat m);
__ZM_EXPORT zCMat zCMatMulNC(const zCMat m1, double k, zCMat m);
__ZM_EXPORT zCMat zCMatDivNC(const zCMat m1, double k, zCMat m);
__ZM_EXPORT zCMat zCMatCMulNC(const zCMat m1, const zComplex *z, zCMat m);
__ZM_EXPORT zCMat zCMatCDivNC(const zCMat m1, const zComplex *z, zCMat m);

__ZM_EXPORT zCMat zCMatAdd(const zCMat m1, const zCMat m2, zCMat m);
__ZM_EXPORT zCMat zCMatSub(const zCMat m1, const zCMat m2, zCMat m);
__ZM_EXPORT zCMat zCMatRev(const zCMat m1, zCMat m);
__ZM_EXPORT zCMat zCMatMul(const zCMat m1, double k, zCMat m);
__ZM_EXPORT zCMat zCMatDiv(const zCMat m1, double k, zCMat m);
__ZM_EXPORT zCMat zCMatCMul(const zCMat m1, const zComplex *z, zCMat m);
__ZM_EXPORT zCMat zCMatCDiv(const zCMat m1, const zComplex *z, zCMat m);

#define zCMatAddNCDRC(m1,m2) zCMatAddNC( (m1), (m2), (m1) )
#define zCMatSubNCDRC(m1,m2) zCMatSubNC( (m1), (m2), (m1) )
#define zCMatRevNCDRC(m)     zCMatRevNC( (m), (m) )
#define zCMatMulNCDRC(m,k)   zCMatMulNC( (m), (k) , (m) )
#define zCMatDivNCDRC(m,k)   zCMatDivNC( (m), (k) , (m) )
#define zCMatCMulNCDRC(m,z)  zCMatCMulNC( (m), (z) , (m) )
#define zCMatCDivNCDRC(m,z)  zCMatCDivNC( (m), (z) , (m) )

#define zCMatAddDRC(m1,m2)   zCMatAdd( (m1), (m2), (m1) )
#define zCMatSubDRC(m1,m2)   zCMatSub( (m1), (m2), (m1) )
#define zCMatRevDRC(m)       zCMatRev( (m), (m) )
#define zCMatMulDRC(m,k)     zCMatMul( (m), (k), (m) )
#define zCMatDivDRC(m,k)     zCMatDiv( (m), (k), (m) )
#define zCMatCMulDRC(m,z)    zCMatMul( (m), (z), (m) )
#define zCMatCDivDRC(m,z)    zCMatDiv( (m), (z), (m) )

/*! \brief multiply a complex matrix and a complex vector, or of two matrices.
 *
 * zCMulMatVecNC() and zCMulMatVec() multiply a column complex vector \a v by a complex matrix \a m,
 * and put the result into \a v.
 * zCMulMatVecNC() does calculation without checking size consistency of vectors and matrices.
 * \return
 * These functions return a pointer to the result.
 */
__ZM_EXPORT zCVec zCMulMatVecNC(const zCMat m, const zCVec v1, zCVec v);
__ZM_EXPORT zCVec zCMulMatVec(const zCMat m, const zCVec v1, zCVec v);

/*! \brief print a complex matrix to a file.
 *
 * zCMatFPrint() prints contents of a complex matrix \a m to current position of the file \a f in the
 * following format:
 *  (r, c) {
 *   z11 z12 ... z1c
 *   z21 z22 ... z2c
 *    .   .  .    .
 *    .   .   .   .
 *    .   .    .  .
 *   zr1 zr2 ... zrc
 *  }
 * where \a r and \a c are row and column size of matrix respectively.
 *
 * zCMatPrint() prints the contents of a complex matrix \a m to the standard output.
 * \return
 * zCMatFPrint() and zCMatPrint() return no values.
 */
__ZM_EXPORT void zCMatFPrint(FILE *fp, const zCMat m);
#define zCMatPrint(m) zCMatFPrint( stdout, (m) )

__END_DECLS

#endif /* __ZM_CMAT_H__ */
