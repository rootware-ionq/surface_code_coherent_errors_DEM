"""Single-qubit gate matrices and identity operators.

Provides the standard Pauli matrices (X, Y, Z) and identity matrices
used as building blocks throughout the simulation.
"""

import numpy as np


def pauli_x() -> np.ndarray:
    """Returns the 2x2 Pauli X matrix.

    X = [[0, 1],
         [1, 0]]
    """
    return np.array([[0, 1], [1, 0]], dtype=complex)


def pauli_y() -> np.ndarray:
    """Returns the 2x2 Pauli Y matrix.

    Y = [[ 0, -i],
         [ i,  0]]
    """
    return np.array([[0, -1j], [1j, 0]], dtype=complex)


def pauli_z() -> np.ndarray:
    """Returns the 2x2 Pauli Z matrix.

    Z = [[1,  0],
         [0, -1]]
    """
    return np.array([[1, 0], [0, -1]], dtype=complex)


def identity_2() -> np.ndarray:
    """Returns the 2x2 identity matrix."""
    return np.eye(2, dtype=complex)


def identity_n(n: int) -> np.ndarray:
    """Returns the 2^n x 2^n identity matrix."""
    return np.eye(1 << n, dtype=complex)
