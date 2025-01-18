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

#define zCMatSizeIsEqual(m1,m2) \
  ( zCMatRowSizeNC(m1) == zCMatRowSizeNC(m2) && \
    zCMatColSizeNC(m1) == zCMatColSizeNC(m2) )

#define zCMatBufNC(m) zArrayBuf(m)
#define zCMatBuf(m)   ( (m) ? zCMatBufNC(m) : NULL )

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

__ZM_EXPORT zCMat zMatToCMat(const zMat m, zCMat cm);

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
 * zCMatMulNC() and zCMatMul() multiply \a m1 by a scalar value \a k, and put the result into \a m.
 *
 * zCMatDivNC() and zCMatDiv() divide \a m1 by \a k, and put the result into \a m.
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
 * \return
 * These functions return a pointer to the result.
 * \notes
 * NC-type functions do calculation without checking the size consistency. If it is not urgent and
 * you are not hasty, you should not use them.
 */
__ZM_EXPORT zCMat zCMatAddNC(const zCMat m1, const zCMat m2, zCMat m);
__ZM_EXPORT zCMat zCMatSubNC(const zCMat m1, const zCMat m2, zCMat m);
__ZM_EXPORT zCMat zCMatRevNC(const zCMat m1, zCMat m);
__ZM_EXPORT zCMat zCMatMulNC(const zCMat m1, const zComplex *z, zCMat m);
__ZM_EXPORT zCMat zCMatDivNC(const zCMat m1, const zComplex *z, zCMat m);

__ZM_EXPORT zCMat zCMatAdd(const zCMat m1, const zCMat m2, zCMat m);
__ZM_EXPORT zCMat zCMatSub(const zCMat m1, const zCMat m2, zCMat m);
__ZM_EXPORT zCMat zCMatRev(const zCMat m1, zCMat m);
__ZM_EXPORT zCMat zCMatMul(const zCMat m1, const zComplex *z, zCMat m);
__ZM_EXPORT zCMat zCMatDiv(const zCMat m1, const zComplex *z, zCMat m);

#define zCMatAddNCDRC(m1,m2) zCMatAddNC( (m1), (m2), (m1) )
#define zCMatSubNCDRC(m1,m2) zCMatSubNC( (m1), (m2), (m1) )
#define zCMatRevNCDRC(m)     zCMatRevNC( (m), (m) )
#define zCMatMulNCDRC(m,z)   zCMatMulNC( (m), (z) , (m) )
#define zCMatDivNCDRC(m,z)   zCMatDivNC( (m), (z) , (m) )

#define zCMatAddDRC(m1,m2)   zCMatAdd( (m1), (m2), (m1) )
#define zCMatSubDRC(m1,m2)   zCMatSub( (m1), (m2), (m1) )
#define zCMatRevDRC(m)       zCMatRev( (m), (m) )
#define zCMatMulDRC(m,z)     zCMatMul( (m), (z), (m) )
#define zCMatDivDRC(m,z)     zCMatDiv( (m), (z), (m) )

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
