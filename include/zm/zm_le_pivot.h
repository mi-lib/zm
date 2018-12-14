/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_pivot - linear equation: pivoting.
 */

#ifndef __ZM_LE_PIVOT_H__
#define __ZM_LE_PIVOT_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* METHOD:
 * zPivoting, zPivotingDiag, zSweepOutVec, zSweepOutMat
 * - matrix pivoting.
 *
 * 'zPivoting()' does partial pivoting on the matrix 'm'.
 * Pivoting is done on 'c'th column.
 * 'index' is a row-ordering index vector. Pivoting is
 * begun from 'r' of 'index' row.
 * For example, suppose 'index' is { 1, 3, 0, 2 } and 'r',
 * 'c' are 2, 1 respectively. This means that the first
 * column of the given matrix has already been pivoted.
 * Then, this function begins pivoting on the second column,
 * examining 0'th and 2'th factor in this order.
 * #
 * 'zPivotingDiag()' only tries to do pivoting on the diagonal
 * factors. Namely, pivoting is done on 'i'th column, begun
 * from 'i' of 'index' row.
 * #
 * 'zSweepOutMat()' ('zSweepOutVec()') sweeps out the
 * factors on 'c'th column of the matrix 'm1' ('m1') other
 * than at 'r'th row and of the vector other than 'r'th,
 * setting them all for 0, and according to it, modifies
 * a matrix 'm2' (a vector 'v').
 * Simultaneously, it divides the factors on 'r'th row of
 * 'm1' ('m').
 * [RETURN VALUE]
 * 'zPivoting()' and 'zPivotingDiag()' returns the index,
 * or the row, of the new pivot on 'c' th column of 'm'.
 * #
 * 'zSweepOutMat()' and 'zSweepOutVec()' return the value
 * originally at the pivot. After calling this function,
 * the value at the same position changes to 1.
 */
__EXPORT int zPivoting(zMat m, zIndex index, int r, int c);
__EXPORT int zPivotingDiag(zMat m, zIndex index, int i);
__EXPORT double zSweepOutVec(zMat m, zVec v, int r, int c);
__EXPORT double zSweepOutMat(zMat m1, zMat m2, int r, int c);

__END_DECLS

#endif /* __ZM_LE_PIVOT_H__ */
