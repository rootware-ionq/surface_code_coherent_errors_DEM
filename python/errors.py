"""Coherent error unitaries for multi-qubit systems.

A coherent Z-rotation error on a single qubit is defined as:

  U(theta) = cos(theta)*I - i*sin(theta)*Z = diag(e^{-i*theta}, e^{+i*theta})

For n qubits the error acts simultaneously on all qubits:

  U_total(theta) = U(theta)^{kron n}

U(theta) is unitary (U*U^dag = I) but not Hermitian for general theta.
"""

import numpy as np


def get_coherent_error_matrix(theta: float, n_qubits: int) -> np.ndarray:
    """Returns the n-qubit coherent error matrix U(theta)^{kron n}.

    The single-qubit factor is U(theta) = cos(theta)*I - i*sin(theta)*Z,
    which equals diag(e^{-i*theta}, e^{+i*theta}).

    The returned matrix is 2^n x 2^n and unitary.
    """
    u = np.array(
        [
            [np.cos(theta) - 1j * np.sin(theta), 0],
            [0, np.cos(theta) + 1j * np.sin(theta)],
        ],
        dtype=complex,
    )
    big_u = u.copy()
    for _ in range(1, n_qubits):
        big_u = np.kron(big_u, u)
    return big_u


def apply_coherent_error(psi: np.ndarray, theta: float) -> np.ndarray:
    """Applies the coherent error U(theta)^{kron n} to a state vector psi.

    Infers the number of qubits from psi.shape[0], which must be a power of 2.
    Returns U(theta) @ |psi>.
    """
    n_qubits = int(np.log2(psi.shape[0]))
    big_u = get_coherent_error_matrix(theta, n_qubits)
    return big_u @ psi
