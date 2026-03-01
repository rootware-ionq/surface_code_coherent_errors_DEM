//! Ancilla measurement probability calculations for the d=3 rotated surface code.
//!
//! Implements syndrome extraction by computing the probability of each possible
//! ancilla measurement outcome for the X and Z stabilizers.

use std::fs::File;
use std::io::{BufWriter, Write};

use faer::{Mat, c64};

use crate::stabilizers::{anti_projection_operator, build_x_stabilizers, projection_operator};

/// Computes the probability of each X ancilla measurement outcome on a given state,
/// and writes the results to `x_ancilla_probs.csv`.
///
/// There are 4 X stabilizers, so there are 2⁴ = 16 possible ancilla patterns
/// A = (a₀, a₁, a₂, a₃) where aᵢ ∈ {0, 1}:
/// - aᵢ = 0 → measure +1 eigenvalue of the i-th X stabilizer (apply projector Pᵢ)
/// - aᵢ = 1 → measure -1 eigenvalue of the i-th X stabilizer (apply anti-projector Qᵢ)
///
/// For each pattern the measurement operator is O = O₀·O₁·O₂·O₃ and the
/// probability is p_A = ⟨φ|O|φ⟩.
///
/// Output file format (CSV, no header):
/// ```text
/// a0,a1,a2,a3,probability
/// ```
/// Readable in Python via `np.loadtxt("x_ancilla_probs.csv", delimiter=",")`.
pub fn perform_x_ancilla_checks(state: &Mat<c64>, filename: Option<&str>) {
    let x_stabs = build_x_stabilizers();

    // Pre-build all projectors and anti-projectors for the 4 X stabilizers.
    let projectors: Vec<(Mat<c64>, Mat<c64>)> = x_stabs
        .iter()
        .map(|s| (projection_operator(s), anti_projection_operator(s)))
        .collect();

    let n_ancilla = projectors.len(); // 4
    let n_patterns = 1usize << n_ancilla; // 16

    let mut results: Vec<([u8; 4], f64)> = Vec::with_capacity(n_patterns);

    for pattern in 0..n_patterns {
        // Decode the bit pattern into (a0, a1, a2, a3).
        let ancilla: [u8; 4] = std::array::from_fn(|i| ((pattern >> i) & 1) as u8);

        // Apply O = O0·O1·O2·O3 sequentially to |φ⟩.
        let mut w = state.clone();
        for (i, (p, q)) in projectors.iter().enumerate() {
            let op = if ancilla[i] == 0 {
                p.as_ref()
            } else {
                q.as_ref()
            };
            w = op * w.as_ref();
        }

        // p_A = ⟨φ|O|φ⟩ = Re(Σᵢ conj(φᵢ) · wᵢ)
        let prob: f64 = (0..state.nrows())
            .map(|i| {
                let phi_i = state.as_ref().get(i, 0);
                let w_i = w.as_ref().get(i, 0);
                phi_i.re * w_i.re + phi_i.im * w_i.im
            })
            .sum();

        results.push((ancilla, prob));
    }

    // Write CSV: a0,a1,a2,a3,probability
    let file;
    if let Some(file_path) = filename {
        file = File::create(file_path).expect("could not create x_ancilla_probs.csv");
    } else {
        file = File::create("x_ancilla_probs.csv").expect("could not create x_ancilla_probs.csv");
    }
    let mut writer = BufWriter::new(file);
    for (ancilla, prob) in &results {
        writeln!(
            writer,
            "{},{},{},{},{:.15e}",
            ancilla[0], ancilla[1], ancilla[2], ancilla[3], prob
        )
        .expect("write failed");
    }

    // Also print to stdout for immediate inspection.
    println!("X ancilla measurement probabilities (a0,a1,a2,a3 -> p):");
    for (ancilla, prob) in &results {
        println!(
            "  ({},{},{},{}) -> {:.6e}",
            ancilla[0], ancilla[1], ancilla[2], ancilla[3], prob
        );
    }
}
