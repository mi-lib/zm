/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lyapnov - linear equation: Lyapnov equation.
 */

#ifndef __ZM_LE_LYAPNOV_H__
#define __ZM_LE_LYAPNOV_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief Lyapnov equation solver.
 *
 * zLELyapnovSolve() solves Lyapnov equation which forms:
 *  X A + A^T X = B
 * where A, B and X are square matrices of the same size.
 * \a a, \a b and \a ans correspond to A, B and X, respectively.
 * \return
 * zLELyapnovSolve() returns the pointer to \a ans if succeeds.
 * The null pointer is returned in the following cases.
 *  1. Any of \a a, \a b and \a ans is not a square matrix.
 *  2. The size of \a a, \a b and \a ans are different from the others.
 *  3. Fails to allocate working memory.
 *  4. \a a is a singular matrix.
 * \notes
 * It internally calls zLESolveGaussDST().
 * \sa
 * zLESolveGaussDST
 */
__ZM_EXPORT zMat zLELyapnovSolve(const zMat a, const zMat b, zMat ans);

__END_DECLS

#endif /* __ZM_LE_LYAPNOV_H__ */
