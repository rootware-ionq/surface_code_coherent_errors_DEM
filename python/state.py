"""Quantum state initialization, logical state preparation, and normalization.

Provides functions to build the initial product state |+>^n, project it
into the d=3 rotated surface code space to obtain the logical |+>_L state,
and renormalize state vectors.
"""

import numpy as np
from stabilizers import build_x_stabilizers, build_z_stabilizers, projection_operator


def plus_ket() -> np.ndarray:
    """Returns the single-qubit |+> = (1/sqrt(2))(|0> + |1>) state as a 2x1 column vector."""
    v = 1.0 / np.sqrt(2.0)
    return np.array([[v], [v]], dtype=complex)


def get_initial_state(n_qubits: int) -> np.ndarray:
    """Returns the n-qubit initial state |+>^n as a 2^n x 1 column vector.

    All qubits start in the |+> state; the full state is formed via successive
    Kronecker products.
    """
    ket = plus_ket()
    state = plus_ket()
    for _ in range(1, n_qubits):
        state = np.kron(state, ket)
    return state


def get_logical_plus_state() -> np.ndarray:
    """Returns the logical |+>_L state of the d=3 rotated surface code (unnormalized).

    Computed by applying all 8 stabilizer projection operators sequentially to
    the initial 9-qubit |+>^9 state:

      |+>_L = (prod_{S in stabilizers} (I + S)/2) |psi_0>

    The result lives in the +1 eigenspace of all stabilizers. Use
    renormalize_state() to obtain a unit-norm state.
    """
    state = get_initial_state(9)
    for s in build_x_stabilizers() + build_z_stabilizers():
        p = projection_operator(s)
        state = p @ state
    return state


def perform_stabilizer(state: np.ndarray) -> np.ndarray:
    """Applies all stabilizer projectors to an arbitrary input state.

    This is the same projection used in get_logical_plus_state() but accepts
    any input state rather than starting from |+>^9.
    """
    state = state.copy()
    for s in build_x_stabilizers() + build_z_stabilizers():
        p = projection_operator(s)
        state = p @ state
    return state


def renormalize_state(state: np.ndarray) -> np.ndarray:
    """Renormalizes a state vector, returning a new unit-norm vector.

    Computes N = sum_n |c_n|^2, then returns state / sqrt(N).
    """
    norm = np.sqrt(np.sum(np.abs(state) ** 2))
    return state / norm
