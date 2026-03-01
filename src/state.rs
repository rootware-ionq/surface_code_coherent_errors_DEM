//! Quantum state initialization, logical state preparation, and normalization.
//!
//! Provides functions to build the initial product state |+⟩^⊗n, project it
//! into the d=3 rotated surface code space to obtain the logical |+⟩_L state,
//! and renormalize state vectors.

use crate::stabilizers::{
    anti_projection_operator, build_x_stabilizers, build_z_stabilizers, projection_operator,
};
use faer::{Mat, c64};

/// Returns the single-qubit |+⟩ = (1/√2)(|0⟩ + |1⟩) state as a 2×1 column vector.
pub fn plus_ket() -> Mat<c64> {
    let v = 1.0_f64 / 2.0_f64.sqrt();
    Mat::from_fn(2, 1, |_, _| c64::new(v, 0.0))
}

/// Returns the n-qubit initial state |+⟩^⊗n as a 2ⁿ×1 column vector.
///
/// All qubits start in the |+⟩ state; the full state is formed via successive
/// Kronecker products.
pub fn get_initial_state(n_qubits: usize) -> Mat<c64> {
    let ket = plus_ket();
    let mut state = plus_ket();
    for _ in 1..n_qubits {
        state = state.as_ref().kron(ket.as_ref());
    }
    state
}

/// Returns the logical |+⟩_L state of the d=3 rotated surface code (unnormalized).
///
/// Computed by applying all 8 stabilizer projection operators sequentially to
/// the initial 9-qubit |+⟩^⊗9 state:
///
/// |+⟩_L = (∏_{S ∈ stabilizers} (I + S)/2) |ψ₀⟩
///
/// The result lives in the +1 eigenspace of all stabilizers. Use
/// [`renormalize_state`] to obtain a unit-norm state.
pub fn get_logical_plus_state() -> Mat<c64> {
    let mut state = get_initial_state(9);
    for s in build_x_stabilizers()
        .iter()
        .chain(build_z_stabilizers().iter())
    {
        let p = projection_operator(s);
        state = p.as_ref() * state.as_ref();
    }
    state
}

/// Applies all stabilizer projectors to an arbitrary input state.
///
/// This is the same projection used in [`get_logical_plus_state`] but accepts
/// any input state rather than starting from |+⟩^⊗9.
pub fn perform_stabilizer(input: &mut Mat<c64>) -> Mat<c64> {
    let mut state = input.clone();
    for s in build_x_stabilizers()
        .iter()
        .chain(build_z_stabilizers().iter())
    {
        let p = projection_operator(s);
        state = p.as_ref() * state.as_ref();
    }
    state
}

/// Renormalizes a state vector, returning a new unit-norm vector.
///
/// Computes 𝒩 = Σₙ |cₙ|², then returns state / √𝒩.
pub fn renormalize_state(state: &mut Mat<c64>) -> Mat<c64> {
    let norm_sq: f64 = (0..state.nrows())
        .map(|i| {
            let v = state.as_ref().get(i, 0);
            v.re * v.re + v.im * v.im
        })
        .sum();
    let norm = norm_sq.sqrt();
    Mat::from_fn(state.nrows(), 1, |i, _| {
        *state.as_ref().get(i, 0) * c64::new(1.0 / norm, 0.0)
    })
}

/// Returns the a modifed initial logical state of the d=3 rotated surface code (unnormalized).
///
/// Computed by applying all 8 stabilizer projection operators sequentially to
/// the initial 9-qubit |+⟩^⊗9 state, except for the Z_2Z_5 stabilizer for which we apply
/// the anti-projection operator
///
/// |psi⟩_L = (I - Z_2Z-5)/2 (∏_{S ∈ stabilizers'} (I + S)/2) |ψ₀⟩
///
/// The result lives in the +1 eigenspace of all stabilizers. Use
/// [`renormalize_state`] to obtain a unit-norm state.
pub fn get_modified_logical_state() -> Mat<c64> {
    let mut state = get_initial_state(9);
    let mut idx = 0;
    for s in build_x_stabilizers()
        .iter()
        .chain(build_z_stabilizers().iter())
    {
        idx += 1;
        let p;
        if idx == 6 {
            p = anti_projection_operator(s);
        } else {
            p = projection_operator(s);
        }
        state = p.as_ref() * state.as_ref();
    }
    state
}
