"""Stabilizer operators and projectors for the d=3 rotated surface code.

The d=3 rotated surface code has 9 data qubits (labelled 0-8) and 8
stabilizers: 4 X-type and 4 Z-type. This module constructs their full
512x512 matrix representations via Kronecker products, and builds the
corresponding projection operators onto the +1 and -1 eigenspaces.
"""

import numpy as np
from gates import identity_2, pauli_x, pauli_z


def x_stabilizer_qubits() -> list[list[int]]:
    """Returns the qubit index sets for each X stabilizer in the d=3 rotated surface code.

    Each inner list gives the qubits on which Pauli X acts; all others get identity.
    The four stabilizers are:
      - X on qubits {0, 1}
      - X on qubits {1, 2, 4, 5}
      - X on qubits {3, 4, 6, 7}
      - X on qubits {7, 8}
    """
    return [[0, 1], [1, 2, 4, 5], [3, 4, 6, 7], [7, 8]]


def z_stabilizer_qubits() -> list[list[int]]:
    """Returns the qubit index sets for each Z stabilizer in the d=3 rotated surface code.

    Each inner list gives the qubits on which Pauli Z acts; all others get identity.
    The four stabilizers are:
      - Z on qubits {0, 1, 3, 4}
      - Z on qubits {2, 5}
      - Z on qubits {3, 6}
      - Z on qubits {4, 5, 7, 8}
    """
    return [[0, 1, 3, 4], [2, 5], [3, 6], [4, 5, 7, 8]]


def build_stabilizer_matrices(
    n_qubits: int,
    stab_qubits: list[list[int]],
    pauli: np.ndarray,
) -> list[np.ndarray]:
    """Builds the full 2^n x 2^n matrix representations of a set of stabilizers.

    For each stabilizer defined by a list of active qubit indices, places `pauli`
    on every active qubit and the 2x2 identity on all others, then Kronecker-products
    across all n_qubits positions to form the full operator.
    """
    id2 = identity_2()
    result = []
    for active in stab_qubits:
        mat = pauli.copy() if 0 in active else id2.copy()
        for qubit in range(1, n_qubits):
            op = pauli if qubit in active else id2
            mat = np.kron(mat, op)
        result.append(mat)
    return result


def build_x_stabilizers() -> list[np.ndarray]:
    """Returns the 4 X-stabilizer matrices for the d=3 rotated surface code.

    Each matrix is 512x512 (2^9 x 2^9). Each stabilizer satisfies S^2 = I.
    """
    return build_stabilizer_matrices(9, x_stabilizer_qubits(), pauli_x())


def build_z_stabilizers() -> list[np.ndarray]:
    """Returns the 4 Z-stabilizer matrices for the d=3 rotated surface code.

    Each matrix is 512x512 (2^9 x 2^9). Each stabilizer satisfies S^2 = I.
    """
    return build_stabilizer_matrices(9, z_stabilizer_qubits(), pauli_z())


def projection_operator(s: np.ndarray) -> np.ndarray:
    """Returns the projection operator P = (I + S) / 2 for a given stabilizer S.

    P projects onto the +1 eigenspace of S. It satisfies P^2 = P (idempotent)
    and P^dag = P (Hermitian), and has the same dimensions as S.
    """
    n = s.shape[0]
    return (np.eye(n, dtype=complex) + s) * 0.5


def anti_projection_operator(s: np.ndarray) -> np.ndarray:
    """Returns the anti-projection operator Q = (I - S) / 2 for a given stabilizer S.

    Q projects onto the -1 eigenspace of S. It satisfies Q^2 = Q (idempotent)
    and Q^dag = Q (Hermitian), and has the same dimensions as S.
    """
    n = s.shape[0]
    return (np.eye(n, dtype=complex) - s) * 0.5
