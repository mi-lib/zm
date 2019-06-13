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

typedef struct{
  int row, col;
  zComplex *elem;
} zCMatStruct;
typedef zCMatStruct * zCMat;

#define _zCMatRowSize(m)         (m)->row
#define _zCMatColSize(m)         (m)->col
#define zCMatRowSize(m)          ( (m) ? _zCMatRowSize(m) : 0 )
#define zCMatColSize(m)          ( (m) ? _zCMatColSize(m) : 0 )
#define zCMatSetRowSize(m,r)     ( _zCMatRowSize(m) = (r) )
#define zCMatSetColSize(m,c)     ( _zCMatColSize(m) = (c) )

#define zCMatSetSize(m,r,c) do{\
  zCMatSetRowSize(m,r);\
  zCMatSetColSize(m,c);\
} while(0)

#define zCMatSizeIsEqual(m1,m2) \
  ( _zCMatRowSize(m1) == _zCMatRowSize(m2) && \
    _zCMatColSize(m1) == _zCMatColSize(m2) )

/*! \brief retrieve and set an element of a complex matrix.
 *
 * zCMatElem() returns the component at \a r th row and
 * \a c th column of a complex matrix \a m.
 *
 * zCMatSetElem() sets the component at \a r th row and
 * \a c th column of a complex matrix \a m for a scalar
 * value \a value.
 * \return
 * zCMatElem() and zCMatSetElem() return the address of
 * the component at \a r th row and \a c th column of \a m.
 */
#define zCMatElem(m,r,c)      ( &(m)->elem[(r)*_zCMatColSize(m)+(c)] )
#define zCMatSetElem(m,r,c,e) zComplexCopy( (e), zCMatElem(m,r,c) )

/*! \brief allocate, free and zero a complex matrix.
 *
 * zCMatAlloc() allocates memory for a new matrix with the size
 * \a row times \a col.
 *
 * zCMatAllocSqr() allocates memory for a new square complex
 * matrix with the size \a size times \a size.
 *
 * zCMatFree() frees the memory allocated for a complex matrix
 * \a m.
 *
 * zCMatClear() clears a complex matrix \a m, namely, set all
 * of components of the matrix for zero.
 * \return
 * zCMatAlloc() and zCMatAllocSqr() return a pointer to the
 * newly allocated memory area.
 *
 * zCMatFree() returns no value.
 *
 * zCMatClear() returns a pointer \a m.
 */
__EXPORT zCMat zCMatAlloc(int row, int col);
#define zCMatAllocSqr(s) zCMatAlloc( (s), (s) )
__EXPORT void zCMatFree(zCMat m);
__EXPORT zCMat zCMatClear(zCMat m);

/*! \brief copy a complex  matrix.
 *
 * zCMatCopyNC() copies a complex matrix \a src to the other
 * \a dest without checking the size consistency between \a src
 * and \a dest.
 *
 * zCMatCopy() copies a complex matrix \a src to the other
 * \a dest.
 *
 * zCMatClone() creates a clone of a complex matrix \a m,
 * namely, allocates memory and copies \a m to it.
 * \return
 * zCMatCopyNC() returns a pointer \a dest.
 *
 * zCMatCopy() returns a pointer \a dest if succeedings, or
 * the null pointer if the size of \a src and \a dest are
 * different.
 *
 * zCMatClone() returns a pointer to the newly created matrix.
 * \notes
 * Since zCMatCopyNC() does not check the size consistency,
 * anything might happen if the size of \a src and \a dest are
 * different. If it is not urgent and you are not hasty, you
 * should use zCMatCopy() for safety.
 */
__EXPORT zCMat zCMatCopyNC(zCMat src, zCMat dest);
__EXPORT zCMat zCMatCopy(zCMat src, zCMat dest);
__EXPORT zCMat zCMatClone(zCMat src);

__EXPORT zCMat zMat2CMat(zMat m, zCMat cm);

/*! \brief check if a complex matrix is tiny.
 *
 * zCMatIsTol() returns the true value if all components of
 * a complex matrix \a m are less than \a tol. Otherwise,
 * the false value is returned.
 *
 * zCMatIsTiny() is the same with zCMatIsTol() except that
 * it compares each component with zTOL, which is defined in
 * "zm_misc.h" as the tolerance.
 * \return
 * zCMatIsTol() and zCMatIsTiny() return a boolean value.
 */
__EXPORT bool zCMatIsTol(zCMat m, double tol);
#define zCMatIsTiny(m) zCMatIsTol( (m), zTOL )

/*! \brief basic arithmetics for the complex matrix.
 *
 * zCMatAddNC() and zCMatAdd() add two complex matrices
 * \a m1 and \a m2, and put the result into \a m.
 *
 * zCMatSubNC() and zCMatSub() subtract \a m2 from \a m1,
 * and put the result into \a m.
 *
 * zCMatRevNC() and zCMatRev() reverse \a m1, and put the
 * result into \a m.
 *
 * zCMatMulNC() and zCMatMul() multiply \a m1 by a scalar
 * value \a k, and put the result into \a m.
 *
 * zCMatDivNC() and zCMatDiv() divide \a m1 by \a k, and
 * put the result into \a m.
 *
 * zCMatAddNCDRC() and zCMatAddDRC() directly add \a m2
 * to \a m1.
 *
 * zCMatSubNCDRC() and zCMatSubDRC() directly subtract
 * \a m2 from \a m1.
 *
 * zCMatRevNCDRC() and zCMatRevDRC() directly reverse \a m.
 *
 * zCMatMulNCDRC() and zCMatMulDRC() directly multiply
 * \a m by \a k.
 *
 * zCMatDivNCDRC() and zCMatDivDRC() directly divide \a m
 * by \a k.
 * \return
 * These functions return a pointer to the result.
 * \notes
 * NC-type functions do calculation without checking the
 * size consistency. If it is not urgent and you are not
 * hasty, you should not use them.
 */
__EXPORT zCMat zCMatAddNC(zCMat m1, zCMat m2, zCMat m);
__EXPORT zCMat zCMatSubNC(zCMat m1, zCMat m2, zCMat m);
__EXPORT zCMat zCMatRevNC(zCMat m1, zCMat m);
__EXPORT zCMat zCMatMulNC(zCMat m1, zComplex *z, zCMat m);
__EXPORT zCMat zCMatDivNC(zCMat m1, zComplex *z, zCMat m);

__EXPORT zCMat zCMatAdd(zCMat m1, zCMat m2, zCMat m);
__EXPORT zCMat zCMatSub(zCMat m1, zCMat m2, zCMat m);
__EXPORT zCMat zCMatRev(zCMat m1, zCMat m);
__EXPORT zCMat zCMatMul(zCMat m1, zComplex *z, zCMat m);
__EXPORT zCMat zCMatDiv(zCMat m1, zComplex *z, zCMat m);

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

/*! \brief multiply a complex matrix and a complex vector,
 * or of two matrices.
 *
 * zCMulMatVecNC() and zCMulMatVec() multiply a column complex
 * vector \a v' by a complex matrix \a m, and put the result
 * into \a v.
 * zCMulMatVecNC() does calculation without checking size
 * consistency of vectors and matrices.
 * \return
 * These functions return a pointer to the result.
 */
__EXPORT zCVec zCMulMatVecNC(zCMat m, zCVec v1, zCVec v);
__EXPORT zCVec zCMulMatVec(zCMat m, zCVec v1, zCVec v);

/*! \brief print a complex matrix to a file.
 *
 * zCMatFPrint() prints contents of a complex matrix \a m to
 * current position of the file \a f' in the following format.
 *  (r, c) {
 *   z11 z12 ... z1c
 *   z21 z22 ... z2c
 *    .   .  .    .
 *    .   .   .   .
 *    .   .    .  .
 *   zr1 zr2 ... zrc
 *  }
 * where \a r and \a c are row and column size of matrix
 * respectively.
 *
 * zCMatPrint() prints the contents of a complex matrix
 * \a m to the standard output.
 * \return
 * zCMatFPrint() and zCMatPrint() return no values.
 */
__EXPORT void zCMatFPrint(FILE *fp, zCMat m);
#define zCMatPrint(m) zCMatFPrint( stdout, (m) )

__END_DECLS

#endif /* __ZM_CMAT_H__ */
