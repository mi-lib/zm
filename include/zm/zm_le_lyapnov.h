/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lyapnov - linear equation: Lyapnov equation.
 */

#ifndef __ZM_LE_LYAPNOV_H__
#define __ZM_LE_LYAPNOV_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* METHOD:
 * zLyapnovSolve
 * - Lyapnov equation solver.
 * [SYNOPSIS]
 * zMat zLyapnovSolve(zMat a, zMat b, zMat ans);
 * [DESCRIPTION]
 * 'zLyapnovSolve()' solves Lyapnov equation which forms:
 *  X A + A^T X = B
 * where A, B and X are square matrices of the same size.
 * 'a', 'b' and 'ans' correspond to A, B and X, respectively.
 * [RETURN VALUE]
 * 'zLyapnovSolve()' returns the pointer to 'ans' if succeeds,
 * or, the null pointer for the following cases.
 *  1. Any of 'a', 'b' and 'ans' is not a square matrix.
 *  2. The size of 'a', 'b' and 'ans' are different from the others.
 *  3. Fails to allocate working memory.
 *  4. 'a' is a singular matrix.
 * [NOTES]
 * It internally calls 'zGaussianSolveDST()'.
 * [SEE ALSO]
 * zGaussianSolveDST
 */
__EXPORT zMat zLyapnovSolve(zMat a, zMat b, zMat ans);

__END_DECLS

#endif /* __ZM_LE_LYAPNOV_H__ */
