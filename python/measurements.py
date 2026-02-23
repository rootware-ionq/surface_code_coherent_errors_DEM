"""Ancilla measurement probability calculations for the d=3 rotated surface code.

Implements syndrome extraction by computing the probability of each possible
ancilla measurement outcome for the X stabilizers.
"""

import numpy as np
from stabilizers import anti_projection_operator, build_x_stabilizers, projection_operator


def perform_x_ancilla_checks(state: np.ndarray) -> list[tuple[tuple[int, ...], float]]:
    """Computes the probability of each X ancilla measurement outcome on a given state,
    and writes the results to x_ancilla_probs.csv.

    There are 4 X stabilizers, so there are 2^4 = 16 possible ancilla patterns
    A = (a0, a1, a2, a3) where a_i in {0, 1}:
      - a_i = 0 -> measure +1 eigenvalue of the i-th X stabilizer (apply projector P_i)
      - a_i = 1 -> measure -1 eigenvalue of the i-th X stabilizer (apply anti-projector Q_i)

    For each pattern the measurement operator is O = O0*O1*O2*O3 and the
    probability is p_A = <phi|O|phi>.

    Output file format (CSV, no header):
      a0,a1,a2,a3,probability
    Readable in Python via np.loadtxt("x_ancilla_probs.csv", delimiter=",").
    """
    x_stabs = build_x_stabilizers()
    projectors = [(projection_operator(s), anti_projection_operator(s)) for s in x_stabs]

    n_ancilla = len(projectors)   # 4
    n_patterns = 1 << n_ancilla   # 16

    results = []
    for pattern in range(n_patterns):
        # Decode the bit pattern into (a0, a1, a2, a3).
        ancilla = tuple((pattern >> i) & 1 for i in range(n_ancilla))

        # Apply O = O0*O1*O2*O3 sequentially to |phi>.
        w = state.copy()
        for i, (p, q) in enumerate(projectors):
            op = p if ancilla[i] == 0 else q
            w = op @ w

        # p_A = <phi|O|phi> = Re(sum_i conj(phi_i) * w_i)
        prob = float(np.real(np.vdot(state.ravel(), w.ravel())))
        results.append((ancilla, prob))

    # Write CSV: a0,a1,a2,a3,probability
    with open("x_ancilla_probs.csv", "w") as f:
        for ancilla, prob in results:
            f.write(f"{ancilla[0]},{ancilla[1]},{ancilla[2]},{ancilla[3]},{prob:.15e}\n")

    # Print to stdout for immediate inspection.
    print("X ancilla measurement probabilities (a0,a1,a2,a3 -> p):")
    for ancilla, prob in results:
        print(f"  {ancilla} -> {prob:.6e}")

    return results
