//! Coherent error unitaries for multi-qubit systems.
//!
//! A coherent Z-rotation error on a single qubit is defined as:
//!
//!   U(θ) = cos(θ)·I − i·sin(θ)·Z = diag(e^{−iθ}, e^{+iθ})
//!
//! For n qubits the error acts simultaneously on all qubits:
//!
//!   𝒰(θ) = U(θ)^{⊗n}
//!
//! U(θ) is unitary (𝒰𝒰† = I) but not Hermitian for general θ.

use faer::{Mat, c64};

/// Returns the n-qubit coherent error matrix 𝒰(θ) = U(θ)^{⊗n}.
///
/// The single-qubit factor is U(θ) = cos(θ)·I − i·sin(θ)·Z,
/// which equals diag(e^{−iθ}, e^{+iθ}).
///
/// The returned matrix is 2ⁿ × 2ⁿ and unitary.
pub fn get_coherent_error_matrix(theta: f64, n_qubits: usize) -> Mat<c64> {
    let cos_t = c64::new(theta.cos(), 0.0);
    let i_sin_t = c64::new(0.0, theta.sin());
    let u = Mat::from_fn(2, 2, |i, j| match (i, j) {
        (0, 0) => cos_t - i_sin_t, // e^{-i theta}
        (1, 1) => cos_t + i_sin_t, // e^{+i theta}
        _ => c64::new(0.0, 0.0),
    });
    let mut big_u = u.clone();
    for _ in 1..n_qubits {
        big_u = big_u.as_ref().kron(u.as_ref());
    }
    big_u
}

/// Applies the coherent error 𝒰(theta) to a state vector `psi`.
///
/// Infers the number of qubits from `psi.nrows()`, which must be a power of 2.
/// Returns 𝒰(theta) · |ψ⟩.
pub fn apply_coherent_error(psi: Mat<c64>, theta: f64) -> Mat<c64> {
    let n_qubits = psi.nrows().ilog2() as usize;
    let big_u = get_coherent_error_matrix(theta, n_qubits);
    big_u.as_ref() * psi.as_ref()
}

/// Same as `apply_coherent_errors` but takes individual thetas for different qubits
pub fn apply_coherent_errors_from_list(psi: Mat<c64>, theta_list: Vec<f64>) -> Mat<c64> {
    let n_qubits = psi.nrows().ilog2() as usize;
    let big_u = get_coherent_error_matrix_from_list(theta_list, n_qubits);
    big_u.as_ref() * psi.as_ref()
}

/// Same as `get_coherent_error_matrix` but takes list of thetas
pub fn get_coherent_error_matrix_from_list(theta_list: Vec<f64>, n_qubits: usize) -> Mat<c64> {
    let cos_t = c64::new(theta_list[0].cos(), 0.0);
    let i_sin_t = c64::new(0.0, theta_list[0].sin());
    let u = Mat::from_fn(2, 2, |i, j| match (i, j) {
        (0, 0) => cos_t - i_sin_t, // e^{-i theta}
        (1, 1) => cos_t + i_sin_t, // e^{+i theta}
        _ => c64::new(0.0, 0.0),
    });
    let mut big_u = u.clone();
    for idx in 1..n_qubits {
        let cos_t = c64::new(theta_list[idx].cos(), 0.0);
        let i_sin_t = c64::new(0.0, theta_list[idx].sin());
        let u = Mat::from_fn(2, 2, |i, j| match (i, j) {
            (0, 0) => cos_t - i_sin_t, // e^{-i theta}
            (1, 1) => cos_t + i_sin_t, // e^{+i theta}
            _ => c64::new(0.0, 0.0),
        });
        big_u = big_u.as_ref().kron(u.as_ref());
    }
    big_u
}
