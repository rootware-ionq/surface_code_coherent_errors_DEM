//! Single-qubit gate matrices and identity operators.
//!
//! Provides the standard Pauli matrices (X, Y, Z) and identity matrices
//! used as building blocks throughout the simulation.

use faer::{Mat, c64};

/// Returns the 2×2 Pauli X matrix.
///
/// ```text
/// X = [[0, 1],
///      [1, 0]]
/// ```
pub fn pauli_x() -> Mat<c64> {
    Mat::from_fn(2, 2, |i, j| match (i, j) {
        (0, 1) | (1, 0) => c64::new(1.0, 0.0),
        _ => c64::new(0.0, 0.0),
    })
}

/// Returns the 2×2 Pauli Y matrix.
///
/// ```text
/// Y = [[ 0, -i],
///      [ i,  0]]
/// ```
pub fn pauli_y() -> Mat<c64> {
    Mat::from_fn(2, 2, |i, j| match (i, j) {
        (0, 1) => c64::new(0.0, -1.0),
        (1, 0) => c64::new(0.0, 1.0),
        _ => c64::new(0.0, 0.0),
    })
}

/// Returns the 2×2 Pauli Z matrix.
///
/// ```text
/// Z = [[1,  0],
///      [0, -1]]
/// ```
pub fn pauli_z() -> Mat<c64> {
    Mat::from_fn(2, 2, |i, j| match (i, j) {
        (0, 0) => c64::new(1.0, 0.0),
        (1, 1) => c64::new(-1.0, 0.0),
        _ => c64::new(0.0, 0.0),
    })
}

/// Returns the 2×2 identity matrix.
pub fn identity_2() -> Mat<c64> {
    Mat::from_fn(2, 2, |i, j| {
        if i == j { c64::new(1.0, 0.0) } else { c64::new(0.0, 0.0) }
    })
}

/// Returns the 2ⁿ × 2ⁿ identity matrix.
pub fn identity_n(n: usize) -> Mat<c64> {
    let dim = 1 << n;
    Mat::from_fn(dim, dim, |i, j| {
        if i == j { c64::new(1.0, 0.0) } else { c64::new(0.0, 0.0) }
    })
}
