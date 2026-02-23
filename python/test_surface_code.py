"""Pytest equivalents of all Rust tests in tests/tests.rs."""

import numpy as np
import pytest

from errors import apply_coherent_error, get_coherent_error_matrix
from gates import identity_2, identity_n, pauli_x, pauli_y, pauli_z
from measurements import perform_x_ancilla_checks
from stabilizers import (
    anti_projection_operator,
    build_x_stabilizers,
    build_z_stabilizers,
    projection_operator,
    x_stabilizer_qubits,
    z_stabilizer_qubits,
)
from state import get_initial_state, get_logical_plus_state, plus_ket, renormalize_state


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def mat_approx_eq(a: np.ndarray, b: np.ndarray, tol: float = 1e-10) -> bool:
    """Returns True if two arrays are element-wise equal within tol."""
    if a.shape != b.shape:
        return False
    return np.all(np.abs(a - b) <= tol)


# ---------------------------------------------------------------------------
# gates
# ---------------------------------------------------------------------------

def test_pauli_x():
    x = pauli_x()
    assert x.shape == (2, 2)
    assert mat_approx_eq(x @ x, identity_2())
    assert abs(x[0, 1] - 1.0) < 1e-10
    assert abs(x[1, 0] - 1.0) < 1e-10
    assert abs(x[0, 0]) < 1e-10
    assert abs(x[1, 1]) < 1e-10


def test_pauli_y():
    y = pauli_y()
    assert y.shape == (2, 2)
    assert mat_approx_eq(y @ y, identity_2())
    assert abs(y[0, 1].imag - (-1.0)) < 1e-10
    assert abs(y[1, 0].imag - 1.0) < 1e-10


def test_pauli_z():
    z = pauli_z()
    assert z.shape == (2, 2)
    assert mat_approx_eq(z @ z, identity_2())
    assert abs(z[0, 0].real - 1.0) < 1e-10
    assert abs(z[1, 1].real - (-1.0)) < 1e-10


# ---------------------------------------------------------------------------
# stabilizers
# ---------------------------------------------------------------------------

def test_x_stabilizer_qubits():
    stabs = x_stabilizer_qubits()
    assert len(stabs) == 4
    assert stabs[0] == [0, 1]
    assert stabs[1] == [1, 2, 4, 5]
    assert stabs[2] == [3, 4, 6, 7]
    assert stabs[3] == [7, 8]


def test_z_stabilizer_qubits():
    stabs = z_stabilizer_qubits()
    assert len(stabs) == 4
    assert stabs[0] == [0, 1, 3, 4]
    assert stabs[1] == [2, 5]
    assert stabs[2] == [3, 6]
    assert stabs[3] == [4, 5, 7, 8]


def test_build_x_stabilizers():
    stabs = build_x_stabilizers()
    assert len(stabs) == 4
    id512 = identity_n(9)
    for s in stabs:
        assert s.shape == (512, 512)
        assert mat_approx_eq(s @ s, id512), "X stabilizer S^2 != I"


def test_build_z_stabilizers():
    stabs = build_z_stabilizers()
    assert len(stabs) == 4
    id512 = identity_n(9)
    for s in stabs:
        assert s.shape == (512, 512)
        assert mat_approx_eq(s @ s, id512), "Z stabilizer S^2 != I"


def test_projection_operator():
    stabs = build_x_stabilizers()
    s = stabs[0]
    p = projection_operator(s)
    assert p.shape == (512, 512)
    assert mat_approx_eq(p @ p, p), "P^2 != P"
    for i, j in [(0, 1), (1, 2), (3, 7), (10, 20)]:
        assert abs(p[i, j].real - p[j, i].real) < 1e-10, f"P not symmetric at re ({i},{j})"
        assert abs(p[i, j].imag + p[j, i].imag) < 1e-10, f"P not Hermitian at im ({i},{j})"


# ---------------------------------------------------------------------------
# anti_projection_operator
# ---------------------------------------------------------------------------

def test_anti_projection_operator():
    stabs = build_x_stabilizers()
    s = stabs[0]
    p = projection_operator(s)
    q = anti_projection_operator(s)

    assert q.shape == (512, 512)

    # Q^2 = Q (idempotent)
    assert mat_approx_eq(q @ q, q), "Q^2 != Q"

    # Q^dag = Q (Hermitian) - spot-check
    for i, j in [(0, 1), (1, 2), (3, 7), (10, 20)]:
        assert abs(q[i, j].real - q[j, i].real) < 1e-10, f"Q not symmetric at re ({i},{j})"
        assert abs(q[i, j].imag + q[j, i].imag) < 1e-10, f"Q not Hermitian at im ({i},{j})"

    # P + Q = I (completeness)
    assert mat_approx_eq(p + q, identity_n(9)), "P + Q != I"

    # P * Q = 0 (orthogonality)
    zero = np.zeros((512, 512), dtype=complex)
    assert mat_approx_eq(p @ q, zero), "P*Q != 0"


# ---------------------------------------------------------------------------
# state
# ---------------------------------------------------------------------------

def test_plus_ket_shape_and_norm():
    ket = plus_ket()
    assert ket.shape == (2, 1)
    norm_sq = float(np.sum(np.abs(ket) ** 2))
    assert abs(norm_sq - 1.0) < 1e-10


def test_get_initial_state():
    state = get_initial_state(3)
    assert state.shape == (8, 1)
    norm_sq = float(np.sum(np.abs(state) ** 2))
    assert abs(norm_sq - 1.0) < 1e-10


def test_get_logical_plus_state():
    state = get_logical_plus_state()
    assert state.shape == (512, 1)
    norm_sq = float(np.sum(np.abs(state) ** 2))
    assert norm_sq > 1e-10, "logical |+> state is zero"
    for s in build_x_stabilizers() + build_z_stabilizers():
        s_state = s @ state
        assert mat_approx_eq(s_state, state), "stabilizer did not fix logical state"


def test_renormalize_state():
    state = np.array([[1 + 1j], [1j], [1.0]], dtype=complex)
    norm_sq = float(np.sum(np.abs(state) ** 2))
    result = renormalize_state(state)
    inv_sqrt_n = 1.0 / np.sqrt(norm_sq)
    expected = np.array([[1 + 1j], [1j], [1.0]], dtype=complex) * inv_sqrt_n
    assert mat_approx_eq(result, expected), "renormalized values incorrect"
    result_norm_sq = float(np.sum(np.abs(result) ** 2))
    assert abs(result_norm_sq - 1.0) < 1e-10, "result not unit norm"


# ---------------------------------------------------------------------------
# errors
# ---------------------------------------------------------------------------

def test_get_coherent_error_matrix():
    n = 3
    dim = 1 << n
    u = get_coherent_error_matrix(0.1 * np.pi, n)
    assert u.shape == (dim, dim)

    # Unitary: U * U^dag = I
    assert mat_approx_eq(u @ u.conj().T, identity_n(n)), "U is not unitary"

    # theta=0 gives identity
    u_zero = get_coherent_error_matrix(0.0, n)
    assert mat_approx_eq(u_zero, identity_n(n)), "U(0) != I"

    # theta=-pi/2: single-qubit U = diag(i, -i) = i*Z
    u_halfpi = get_coherent_error_matrix(-np.pi / 2, 1)
    assert abs(u_halfpi[0, 0].real) < 1e-10
    assert abs(u_halfpi[0, 0].imag - 1.0) < 1e-10
    assert abs(u_halfpi[1, 1].real) < 1e-10
    assert abs(u_halfpi[1, 1].imag - (-1.0)) < 1e-10
    assert abs(u_halfpi[0, 1].real) < 1e-10
    assert abs(u_halfpi[1, 0].real) < 1e-10


def test_apply_coherent_error():
    logical = get_logical_plus_state()
    psi = renormalize_state(logical)

    errored = apply_coherent_error(psi.copy(), 0.1 * np.pi)
    assert errored.shape == (512, 1)

    # theta=0 leaves state unchanged
    unchanged = apply_coherent_error(psi.copy(), 0.0)
    assert mat_approx_eq(unchanged, psi), "U(0) did not preserve state"

    # Norm is preserved
    norm_sq = float(np.sum(np.abs(errored) ** 2))
    assert abs(norm_sq - 1.0) < 1e-10, "norm not preserved after error"


# ---------------------------------------------------------------------------
# measurements
# ---------------------------------------------------------------------------

def test_perform_x_ancilla_checks():
    logical = get_logical_plus_state()
    psi = renormalize_state(logical)
    state = apply_coherent_error(psi.copy(), 0.1 * np.pi)

    x_stabs = build_x_stabilizers()
    projectors = [(projection_operator(s), anti_projection_operator(s)) for s in x_stabs]

    n_patterns = 1 << len(projectors)
    prob_sum = 0.0

    for pattern in range(n_patterns):
        ancilla = tuple((pattern >> i) & 1 for i in range(len(projectors)))
        w = state.copy()
        for i, (p, q) in enumerate(projectors):
            op = p if ancilla[i] == 0 else q
            w = op @ w
        prob = float(np.real(np.vdot(state.ravel(), w.ravel())))
        assert prob >= -1e-10, f"negative probability at pattern {pattern}: {prob}"
        prob_sum += prob

    assert abs(prob_sum - 1.0) < 1e-10, f"probabilities do not sum to 1: {prob_sum}"

    # Exercise the full function (also writes x_ancilla_probs.csv).
    perform_x_ancilla_checks(state)
