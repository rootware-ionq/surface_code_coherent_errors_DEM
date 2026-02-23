# Surface Code Simulation — Project Summary

## Overview

This project implements a simulation of the **d=3 rotated surface code** under coherent Z-rotation errors, first in Rust (using `faer` for linear algebra) and then ported to Python (using NumPy). The goal is to compute the probability distribution over all possible X ancilla measurement outcomes after applying a coherent error to the logical `|+>_L` state.

---

## Step 1 — Building Blocks (`src/gates.rs`, `src/stabilizers.rs`, `src/state.rs`)

### State Initialization
- **`plus_ket()`**: Returns the single-qubit `|+> = (1/√2)(|0> + |1>)` state as a 2×1 complex column vector.
- **`get_initial_state(n_qubits)`**: Returns the n-qubit product state `|+>^⊗n` as a 2ⁿ×1 vector via successive Kronecker products.

### Pauli Matrices
- **`pauli_x()`**, **`pauli_y()`**, **`pauli_z()`**: Standard 2×2 Pauli matrices as complex matrices.
- **`identity_2()`**, **`identity_n(n)`**: 2×2 and 2ⁿ×2ⁿ identity matrices.

### d=3 Rotated Surface Code Stabilizers
The d=3 rotated surface code has 9 data qubits (labelled 0–8) and 8 stabilizers:

- **X stabilizers** act on qubit sets: `{0,1}`, `{1,2,4,5}`, `{3,4,6,7}`, `{7,8}`
- **Z stabilizers** act on qubit sets: `{0,1,3,4}`, `{2,5}`, `{3,6}`, `{4,5,7,8}`

- **`x_stabilizer_qubits()`** / **`z_stabilizer_qubits()`**: Return these index sets as `Vec<Vec<usize>>`.
- **`build_stabilizer_matrices(n_qubits, stab_qubits, pauli)`**: Generic Kronecker-product builder — places `pauli` on active qubits and `I` on the rest, producing a 2ⁿ×2ⁿ matrix.
- **`build_x_stabilizers()`** / **`build_z_stabilizers()`**: Return the 4 X and 4 Z stabilizer matrices (512×512 each). Each satisfies S² = I.

---

## Step 2 — Logical State Preparation (`src/stabilizers.rs`, `src/state.rs`)

### Projection Operators
- **`projection_operator(S)`**: Returns `P = (I + S) / 2`, projecting onto the +1 eigenspace of stabilizer S. Satisfies P² = P and P† = P.

### Logical |+>_L State
- **`get_logical_plus_state()`**: Prepares the logical codeword state by applying all 8 stabilizer projectors sequentially to `|+>^⊗9`:

  ```
  |+>_L = (∏_{S ∈ stabilizers} (I + S)/2) |ψ₀>
  ```

  Result: 32 non-zero amplitudes of equal magnitude `1/(4√2)`, all real — a uniform superposition over the 32 codeword basis states.

---

## Step 3 — Normalization and Coherent Errors (`src/state.rs`, `src/errors.rs`)

### Normalization
- **`renormalize_state(state)`**: Divides the state by √(Σ|cₙ|²) to produce a unit-norm vector.

### Coherent Error Unitary
The single-qubit coherent error is a Z-rotation:
```
U(θ) = cos(θ)·I − i·sin(θ)·Z = diag(e^{−iθ}, e^{+iθ})
```
Acting on all n qubits simultaneously:
```
𝒰(θ) = U(θ)^{⊗n}
```
`𝒰` is **unitary** (𝒰𝒰† = I) but not Hermitian for general θ.

- **`get_coherent_error_matrix(theta, n_qubits)`**: Constructs the full 2ⁿ×2ⁿ error unitary via Kronecker products.
- **`apply_coherent_error(psi, theta)`**: Infers n_qubits from `psi.nrows()` (must be power of 2), then returns `𝒰(θ)|ψ>`.

### In `main()`
With `θ = 0.1π`, the coherent error preserves all 32 basis-state populations (norm² = 1/32 each) — Z-rotation only adds phases, not amplitude mixing.

---

## Step 4 — Ancilla Measurements (`src/stabilizers.rs`, `src/measurements.rs`)

### Anti-Projection Operator
- **`anti_projection_operator(S)`**: Returns `Q = (I − S) / 2`, projecting onto the **−1** eigenspace of S. Satisfies Q² = Q, Q† = Q, P + Q = I, P·Q = 0.

### X Ancilla Checks
- **`perform_x_ancilla_checks(state)`**: Iterates over all 16 ancilla patterns A = (a₀, a₁, a₂, a₃) where aᵢ ∈ {0, 1}:
  - aᵢ = 0 → apply projector Pᵢ (+1 measurement outcome)
  - aᵢ = 1 → apply anti-projector Qᵢ (−1 measurement outcome)
  - Composite operator: O = O₀·O₁·O₂·O₃ applied to |φ>
  - Probability: p_A = ⟨φ|O|φ⟩ = Re(Σᵢ conj(φᵢ)·(O|φ>)ᵢ)

  Writes results to `x_ancilla_probs.csv` (NumPy-readable: `np.loadtxt(..., delimiter=",")`).

### Results (θ = 0.1π)
| Ancilla (a₀,a₁,a₂,a₃) | Probability |
|------------------------|-------------|
| (0,0,0,0) | 0.3094 |
| (0,1,0,0) | 0.1555 |
| (0,0,1,0) | 0.1555 |
| (0,1,1,0) | 0.0640 |
| (1,0,1,0) | 0.0465 |
| ... | ... |
| (1,1,1,1) | 0.0075 |

The (0,0,0,0) pattern (all +1 outcomes, no error detected) dominates at ~31%, as expected for a small rotation angle.

---

## Step 5 — Python Port (`python/`)

A complete NumPy-based Python implementation mirroring the Rust code:

| File | Contents |
|------|----------|
| `gates.py` | `pauli_x`, `pauli_y`, `pauli_z`, `identity_2`, `identity_n` |
| `stabilizers.py` | Stabilizer indices, matrix builders, `projection_operator`, `anti_projection_operator` |
| `state.py` | `plus_ket`, `get_initial_state`, `get_logical_plus_state`, `perform_stabilizer`, `renormalize_state` |
| `errors.py` | `get_coherent_error_matrix`, `apply_coherent_error` |
| `measurements.py` | `perform_x_ancilla_checks` → writes `x_ancilla_probs.csv` |
| `main.py` | Entry point mirroring `main.rs` |
| `test_surface_code.py` | 16 pytest tests, one per Rust test (all passing) |

**Agreement with Rust**: CSV outputs match to within **1.9×10⁻¹⁶** (machine epsilon), with all 16 probabilities agreeing at 6+ significant figures.

---

## Code Structure

```
test_claude/
├── Cargo.toml                  ← faer = "0.24" dependency
├── src/
│   ├── lib.rs                  ← crate root; declares pub modules
│   ├── main.rs                 ← thin entry point
│   ├── gates.rs                ← Pauli matrices, identity operators
│   ├── stabilizers.rs          ← stabilizer construction, projectors
│   ├── state.rs                ← state prep, normalization
│   ├── errors.rs               ← coherent error unitaries
│   └── measurements.rs         ← ancilla probability calculations
├── tests/
│   └── tests.rs                ← 16 integration tests (all passing)
└── python/
    ├── gates.py
    ├── stabilizers.py
    ├── state.py
    ├── errors.py
    ├── measurements.py
    ├── main.py
    └── test_surface_code.py    ← 16 pytest tests (all passing)
```

---

## Key Physics Notes

1. **Stabilizer code**: The d=3 rotated surface code encodes 1 logical qubit in 9 physical qubits using 8 stabilizer constraints (4 X-type, 4 Z-type).

2. **Logical state preparation**: Starting from `|+>^⊗9` and projecting with all stabilizers onto the +1 eigenspace yields the logical `|+>_L` — a uniform superposition of 32 computational basis states.

3. **Coherent error**: U(θ) = cos(θ)I − i·sin(θ)Z is a Z-rotation. It commutes with Z stabilizers but anti-commutes with X stabilizers, introducing correlations in the X syndrome measurements.

4. **Syndrome extraction**: The X ancilla measurement distribution encodes information about whether errors occurred, enabling error correction. The deviation of probabilities from the error-free (0,0,0,0) = 1 case quantifies the effect of the coherent error.
