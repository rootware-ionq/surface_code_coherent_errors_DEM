"""Entry point: replicates main.rs in Python.

Prepares the logical |+>_L state, applies coherent errors, and runs
X ancilla measurements, writing results to x_ancilla_probs.csv.
"""

import numpy as np
from errors import apply_coherent_error
from measurements import perform_x_ancilla_checks
from state import get_logical_plus_state, renormalize_state


def main():
    logical = get_logical_plus_state()
    logical_plus_state = renormalize_state(logical)

    theta = 0.1 * np.pi
    coherent_error_state = apply_coherent_error(logical_plus_state, theta)

    perform_x_ancilla_checks(coherent_error_state)


if __name__ == "__main__":
    main()
