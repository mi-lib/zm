/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_pivot - linear equation: pivoting.
 */

#ifndef __ZM_LE_PIVOT_H__
#define __ZM_LE_PIVOT_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief matrix pivoting.
 *
 * zMatPivoting() does partial pivoting on the matrix \a m on the \a c th column.
 * \a index is a row-ordering index vector.
 * Pivoting begins from \a r of \a index row. For example, if \a index is { 1, 3, 0, 2 } and \a r and
 * \a c are 2 and 1, respectively, the first column of \a m has to be already pivoted. Then, this
 * function begins pivoting on the second column, examining 0th and 2nd factor in this order.
 *
 * zMatPivotingDiag() only tries to do pivoting on the diagonal components.
 * Namely, pivoting is done on the \a i th column, beginning from \a i of index row.
 *
 * zMatMatSweepOut() and zMatVecSweepOut() sweep out components on the \a c th
 * column of the matrix \a m1 other than at the \a r th row and of the
 * vector other than \a r th, setting them all for 0, and according to
 * it, modifies the matrix \a m2 and the vector \a v, respectively.
 * They simultaneously divide the components on the \a r th row of \a m1
 * and \a m, respectively.
 * \return
 * zMatPivoting() and zMatPivotingDiag() return the index, or the row, of the
 * new pivot on the \a c th column of \a m.
 *
 * zMatMatSweepOut() and zMatVecSweepOut() return the value originally at the
 * pivot. After calling these functions, the value at the same position
 * changes to 1.
 */
__ZM_EXPORT int zMatPivoting(zMat m, zIndex index, int r, int c);
__ZM_EXPORT int zMatPivotingDiag(zMat m, zIndex index, int i);
__ZM_EXPORT double zMatMatSweepOut(zMat m1, zMat m2, int r, int c);
__ZM_EXPORT double zMatVecSweepOut(zMat m, zVec v, int r, int c);

__END_DECLS

#endif /* __ZM_LE_PIVOT_H__ */
