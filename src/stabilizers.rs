//! Stabilizer operators and projectors for the d=3 rotated surface code.
//!
//! The d=3 rotated surface code has 9 data qubits (labelled 0–8) and 8
//! stabilizers: 4 X-type and 4 Z-type. This module constructs their full
//! 512×512 matrix representations via Kronecker products, and builds the
//! corresponding projection operators onto the +1 eigenspace.

use crate::gates::{identity_2, pauli_x, pauli_z};
use faer::{Mat, c64};

/// Returns the qubit index sets for each X stabilizer in the d=3 rotated surface code.
///
/// Each inner `Vec<usize>` lists the qubits on which Pauli X acts; all other
/// qubits receive the identity. The four stabilizers are:
/// - X on qubits {0, 1}
/// - X on qubits {1, 2, 4, 5}
/// - X on qubits {3, 4, 6, 7}
/// - X on qubits {7, 8}
pub fn x_stabilizer_qubits() -> Vec<Vec<usize>> {
    vec![vec![0, 1], vec![1, 2, 4, 5], vec![3, 4, 6, 7], vec![7, 8]]
}

/// Returns the qubit index sets for each Z stabilizer in the d=3 rotated surface code.
///
/// Each inner `Vec<usize>` lists the qubits on which Pauli Z acts; all other
/// qubits receive the identity. The four stabilizers are:
/// - Z on qubits {0, 1, 3, 4}
/// - Z on qubits {2, 5}
/// - Z on qubits {3, 6}
/// - Z on qubits {4, 5, 7, 8}
pub fn z_stabilizer_qubits() -> Vec<Vec<usize>> {
    vec![vec![0, 1, 3, 4], vec![2, 5], vec![3, 6], vec![4, 5, 7, 8]]
}

/// Builds the full 2ⁿ × 2ⁿ matrix representations of a set of stabilizers.
///
/// For each stabilizer defined by a list of active qubit indices, this function
/// places `pauli` on every active qubit and the 2×2 identity on all others,
/// then Kronecker-products across all `n_qubits` positions to form the full operator.
pub fn build_stabilizer_matrices(
    n_qubits: usize,
    stab_qubits: &[Vec<usize>],
    pauli: &Mat<c64>,
) -> Vec<Mat<c64>> {
    let id = identity_2();
    stab_qubits
        .iter()
        .map(|active| {
            let first_op = if active.contains(&0) {
                pauli.as_ref()
            } else {
                id.as_ref()
            };
            let mut mat = first_op.to_owned();
            for qubit in 1..n_qubits {
                let op = if active.contains(&qubit) {
                    pauli.as_ref()
                } else {
                    id.as_ref()
                };
                mat = mat.as_ref().kron(op);
            }
            mat
        })
        .collect()
}

/// Returns the 4 X-stabilizer matrices for the d=3 rotated surface code.
///
/// Each matrix is 512×512 (2⁹ × 2⁹). Each stabilizer satisfies S² = I.
pub fn build_x_stabilizers() -> Vec<Mat<c64>> {
    build_stabilizer_matrices(9, &x_stabilizer_qubits(), &pauli_x())
}

/// Returns the 4 Z-stabilizer matrices for the d=3 rotated surface code.
///
/// Each matrix is 512×512 (2⁹ × 2⁹). Each stabilizer satisfies S² = I.
pub fn build_z_stabilizers() -> Vec<Mat<c64>> {
    build_stabilizer_matrices(9, &z_stabilizer_qubits(), &pauli_z())
}

/// Returns the projection operator P = (I + S) / 2 for a given stabilizer S.
///
/// P projects onto the +1 eigenspace of S. It satisfies P² = P (idempotent)
/// and P† = P (Hermitian), and has the same dimensions as S.
pub fn projection_operator(s: &Mat<c64>) -> Mat<c64> {
    let n = s.nrows();
    Mat::from_fn(n, n, |i, j| {
        let id_val = if i == j {
            c64::new(1.0, 0.0)
        } else {
            c64::new(0.0, 0.0)
        };
        (id_val + *s.as_ref().get(i, j)) * c64::new(0.5, 0.0)
    })
}

/// Returns the anti-projection operator Q = (I - S) / 2 for a given stabilizer S.
///
/// Q projects onto the -1 eigenspace of S. It satisfies Q² = Q (idempotent)
/// and Q† = Q (Hermitian), and has the same dimensions as S.
pub fn anti_projection_operator(s: &Mat<c64>) -> Mat<c64> {
    let n = s.nrows();
    Mat::from_fn(n, n, |i, j| {
        let id_val = if i == j {
            c64::new(1.0, 0.0)
        } else {
            c64::new(0.0, 0.0)
        };
        (id_val - *s.as_ref().get(i, j)) * c64::new(0.5, 0.0)
    })
}
