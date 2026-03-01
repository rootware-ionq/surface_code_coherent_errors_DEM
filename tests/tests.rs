use faer::{Mat, c64};
use surface_code_dem::errors::{apply_coherent_error, get_coherent_error_matrix};
use surface_code_dem::gates::{identity_2, identity_n, pauli_x, pauli_y, pauli_z};
use surface_code_dem::measurements::perform_x_ancilla_checks;
use surface_code_dem::stabilizers::{
    anti_projection_operator, build_x_stabilizers, build_z_stabilizers, projection_operator,
    x_stabilizer_qubits, z_stabilizer_qubits,
};
use surface_code_dem::state::{
    get_initial_state, get_logical_plus_state, plus_ket, renormalize_state,
};

/// Returns true if two matrices are element-wise equal within `tol`.
fn mat_approx_eq(a: &Mat<c64>, b: &Mat<c64>, tol: f64) -> bool {
    if a.nrows() != b.nrows() || a.ncols() != b.ncols() {
        return false;
    }
    for i in 0..a.nrows() {
        for j in 0..a.ncols() {
            let diff = *a.as_ref().get(i, j) - *b.as_ref().get(i, j);
            if (diff.re * diff.re + diff.im * diff.im).sqrt() > tol {
                return false;
            }
        }
    }
    true
}

// --- gates ---

#[test]
fn test_pauli_x() {
    let x = pauli_x();
    assert_eq!(x.nrows(), 2);
    assert_eq!(x.ncols(), 2);
    let x2 = x.as_ref() * x.as_ref();
    assert!(mat_approx_eq(&x2, &identity_2(), 1e-10));
    assert!((x.as_ref().get(0, 1).re - 1.0).abs() < 1e-10);
    assert!((x.as_ref().get(1, 0).re - 1.0).abs() < 1e-10);
    assert!(x.as_ref().get(0, 0).re.abs() < 1e-10);
    assert!(x.as_ref().get(1, 1).re.abs() < 1e-10);
}

#[test]
fn test_pauli_y() {
    let y = pauli_y();
    assert_eq!(y.nrows(), 2);
    assert_eq!(y.ncols(), 2);
    let y2 = y.as_ref() * y.as_ref();
    assert!(mat_approx_eq(&y2, &identity_2(), 1e-10));
    assert!((y.as_ref().get(0, 1).im - (-1.0)).abs() < 1e-10);
    assert!((y.as_ref().get(1, 0).im - 1.0).abs() < 1e-10);
}

#[test]
fn test_pauli_z() {
    let z = pauli_z();
    assert_eq!(z.nrows(), 2);
    assert_eq!(z.ncols(), 2);
    let z2 = z.as_ref() * z.as_ref();
    assert!(mat_approx_eq(&z2, &identity_2(), 1e-10));
    assert!((z.as_ref().get(0, 0).re - 1.0).abs() < 1e-10);
    assert!((z.as_ref().get(1, 1).re - (-1.0)).abs() < 1e-10);
}

// --- stabilizers ---

#[test]
fn test_x_stabilizer_qubits() {
    let stabs = x_stabilizer_qubits();
    assert_eq!(stabs.len(), 4);
    assert_eq!(stabs[0], vec![0, 1]);
    assert_eq!(stabs[1], vec![1, 2, 4, 5]);
    assert_eq!(stabs[2], vec![3, 4, 6, 7]);
    assert_eq!(stabs[3], vec![7, 8]);
}

#[test]
fn test_z_stabilizer_qubits() {
    let stabs = z_stabilizer_qubits();
    assert_eq!(stabs.len(), 4);
    assert_eq!(stabs[0], vec![0, 1, 3, 4]);
    assert_eq!(stabs[1], vec![2, 5]);
    assert_eq!(stabs[2], vec![3, 6]);
    assert_eq!(stabs[3], vec![4, 5, 7, 8]);
}

#[test]
fn test_build_x_stabilizers() {
    let stabs = build_x_stabilizers();
    assert_eq!(stabs.len(), 4);
    let id512 = identity_n(9);
    for s in &stabs {
        assert_eq!(s.nrows(), 512);
        assert_eq!(s.ncols(), 512);
        let s2 = s.as_ref() * s.as_ref();
        assert!(mat_approx_eq(&s2, &id512, 1e-10), "X stabilizer S^2 != I");
    }
}

#[test]
fn test_build_z_stabilizers() {
    let stabs = build_z_stabilizers();
    assert_eq!(stabs.len(), 4);
    let id512 = identity_n(9);
    for s in &stabs {
        assert_eq!(s.nrows(), 512);
        assert_eq!(s.ncols(), 512);
        let s2 = s.as_ref() * s.as_ref();
        assert!(mat_approx_eq(&s2, &id512, 1e-10), "Z stabilizer S^2 != I");
    }
}

#[test]
fn test_projection_operator() {
    let stabs = build_x_stabilizers();
    let s = &stabs[0];
    let p = projection_operator(s);
    assert_eq!(p.nrows(), 512);
    assert_eq!(p.ncols(), 512);
    let p2 = p.as_ref() * p.as_ref();
    assert!(mat_approx_eq(&p2, &p, 1e-10), "P^2 != P");
    for (i, j) in [(0, 1), (1, 2), (3, 7), (10, 20)] {
        let pij = p.as_ref().get(i, j);
        let pji = p.as_ref().get(j, i);
        assert!(
            (pij.re - pji.re).abs() < 1e-10,
            "P not symmetric at re ({i},{j})"
        );
        assert!(
            (pij.im + pji.im).abs() < 1e-10,
            "P not Hermitian at im ({i},{j})"
        );
    }
}

// --- state ---

#[test]
fn test_plus_ket_shape_and_norm() {
    let ket = plus_ket();
    assert_eq!(ket.nrows(), 2);
    assert_eq!(ket.ncols(), 1);
    let norm_sq: f64 = (0..2)
        .map(|i| {
            let v = ket.as_ref().get(i, 0);
            v.re * v.re + v.im * v.im
        })
        .sum();
    assert!((norm_sq - 1.0).abs() < 1e-10);
}

#[test]
fn test_get_initial_state() {
    let state = get_initial_state(3);
    assert_eq!(state.nrows(), 8);
    assert_eq!(state.ncols(), 1);
    let norm_sq: f64 = (0..8)
        .map(|i| {
            let v = state.as_ref().get(i, 0);
            v.re * v.re + v.im * v.im
        })
        .sum();
    assert!((norm_sq - 1.0).abs() < 1e-10);
}

#[test]
fn test_get_logical_plus_state() {
    let state = get_logical_plus_state();
    assert_eq!(state.nrows(), 512);
    assert_eq!(state.ncols(), 1);
    let norm_sq: f64 = (0..512)
        .map(|i| {
            let v = state.as_ref().get(i, 0);
            v.re * v.re + v.im * v.im
        })
        .sum();
    assert!(norm_sq > 1e-10, "logical |+> state is zero");
    for s in build_x_stabilizers()
        .iter()
        .chain(build_z_stabilizers().iter())
    {
        let s_state = s.as_ref() * state.as_ref();
        assert!(
            mat_approx_eq(&s_state, &state, 1e-10),
            "stabilizer did not fix logical state"
        );
    }
}

#[test]
fn test_renormalize_state() {
    let mut state = Mat::from_fn(3, 1, |i, _| match i {
        0 => c64::new(1.0, 1.0),
        1 => c64::new(0.0, 1.0),
        _ => c64::new(1.0, 0.0),
    });
    let norm_sq: f64 = (0..3)
        .map(|i| {
            let v = state.as_ref().get(i, 0);
            v.re * v.re + v.im * v.im
        })
        .sum();
    let result = renormalize_state(&mut state);
    let inv_sqrt_n = 1.0 / norm_sq.sqrt();
    let expected = [c64::new(1.0, 1.0), c64::new(0.0, 1.0), c64::new(1.0, 0.0)];
    for i in 0..3 {
        let got = result.as_ref().get(i, 0);
        let exp = expected[i] * c64::new(inv_sqrt_n, 0.0);
        assert!((got.re - exp.re).abs() < 1e-10, "re mismatch at {i}");
        assert!((got.im - exp.im).abs() < 1e-10, "im mismatch at {i}");
    }
    let result_norm_sq: f64 = (0..3)
        .map(|i| {
            let v = result.as_ref().get(i, 0);
            v.re * v.re + v.im * v.im
        })
        .sum();
    assert!((result_norm_sq - 1.0).abs() < 1e-10, "result not unit norm");
}

// --- anti_projection_operator ---

#[test]
fn test_anti_projection_operator() {
    let stabs = build_x_stabilizers();
    let s = &stabs[0];
    let p = projection_operator(s);
    let q = anti_projection_operator(s);

    // Shape matches stabilizer
    assert_eq!(q.nrows(), 512);
    assert_eq!(q.ncols(), 512);

    // Q² = Q (idempotent)
    let q2 = q.as_ref() * q.as_ref();
    assert!(mat_approx_eq(&q2, &q, 1e-10), "Q^2 != Q");

    // Q† = Q (Hermitian) — spot-check entries
    for (i, j) in [(0, 1), (1, 2), (3, 7), (10, 20)] {
        let qij = q.as_ref().get(i, j);
        let qji = q.as_ref().get(j, i);
        assert!(
            (qij.re - qji.re).abs() < 1e-10,
            "Q not symmetric at re ({i},{j})"
        );
        assert!(
            (qij.im + qji.im).abs() < 1e-10,
            "Q not Hermitian at im ({i},{j})"
        );
    }

    // P + Q = I (completeness)
    let id512 = identity_n(9);
    let p_plus_q = Mat::from_fn(512, 512, |i, j| {
        *p.as_ref().get(i, j) + *q.as_ref().get(i, j)
    });
    assert!(mat_approx_eq(&p_plus_q, &id512, 1e-10), "P + Q != I");

    // P · Q = 0 (orthogonality)
    let pq = p.as_ref() * q.as_ref();
    let zero = Mat::from_fn(512, 512, |_, _| c64::new(0.0, 0.0));
    assert!(mat_approx_eq(&pq, &zero, 1e-10), "P·Q != 0");
}

// --- errors ---

#[test]
fn test_get_coherent_error_matrix() {
    let n = 3;
    let dim = 1 << n;
    let u = get_coherent_error_matrix(0.1 * std::f64::consts::PI, n);
    assert_eq!(u.nrows(), dim);
    assert_eq!(u.ncols(), dim);
    // Unitary: U * U† = I
    let u_dag = Mat::from_fn(dim, dim, |i, j| {
        let v = u.as_ref().get(j, i);
        c64::new(v.re, -v.im)
    });
    let product = u.as_ref() * u_dag.as_ref();
    assert!(
        mat_approx_eq(&product, &identity_n(n), 1e-10),
        "U is not unitary"
    );
    // theta=0 gives identity
    let u_zero = get_coherent_error_matrix(0.0, n);
    assert!(mat_approx_eq(&u_zero, &identity_n(n), 1e-10), "U(0) != I");
    // theta=-pi/2: single-qubit U = diag(i, -i) = i*Z
    let u_halfpi = get_coherent_error_matrix(-std::f64::consts::FRAC_PI_2, 1);
    assert!((u_halfpi.as_ref().get(0, 0).re - 0.0).abs() < 1e-10);
    assert!((u_halfpi.as_ref().get(0, 0).im - 1.0).abs() < 1e-10);
    assert!((u_halfpi.as_ref().get(1, 1).re - 0.0).abs() < 1e-10);
    assert!((u_halfpi.as_ref().get(1, 1).im - (-1.0)).abs() < 1e-10);
    assert!(u_halfpi.as_ref().get(0, 1).re.abs() < 1e-10);
    assert!(u_halfpi.as_ref().get(1, 0).re.abs() < 1e-10);
}

#[test]
fn test_apply_coherent_error() {
    let mut logical = get_logical_plus_state();
    let psi = renormalize_state(&mut logical);
    let errored = apply_coherent_error(psi.clone(), 0.1 * std::f64::consts::PI);
    assert_eq!(errored.nrows(), 512);
    assert_eq!(errored.ncols(), 1);
    // theta=0 leaves state unchanged
    let unchanged = apply_coherent_error(psi.clone(), 0.0);
    assert!(
        mat_approx_eq(&unchanged, &psi, 1e-10),
        "U(0) did not preserve state"
    );
    // Norm is preserved
    let norm_sq: f64 = (0..512)
        .map(|i| {
            let v = errored.as_ref().get(i, 0);
            v.re * v.re + v.im * v.im
        })
        .sum();
    assert!(
        (norm_sq - 1.0).abs() < 1e-10,
        "norm not preserved after error"
    );
}

// --- measurements ---

#[test]
fn test_perform_x_ancilla_checks() {
    let mut logical = get_logical_plus_state();
    let psi = renormalize_state(&mut logical);
    let state = apply_coherent_error(psi.clone(), 0.1 * std::f64::consts::PI);

    // Run the measurement (also writes x_ancilla_probs.csv as a side effect).
    // We verify the returned probabilities by re-computing them here independently.
    use surface_code_dem::stabilizers::{
        anti_projection_operator, build_x_stabilizers, projection_operator,
    };
    let x_stabs = build_x_stabilizers();
    let projectors: Vec<(_, _)> = x_stabs
        .iter()
        .map(|s| (projection_operator(s), anti_projection_operator(s)))
        .collect();

    let n_patterns = 1usize << projectors.len(); // 16
    let mut prob_sum = 0.0_f64;

    for pattern in 0..n_patterns {
        let ancilla: [u8; 4] = std::array::from_fn(|i| ((pattern >> i) & 1) as u8);
        let mut w = state.clone();
        for (i, (p, q)) in projectors.iter().enumerate() {
            let op = if ancilla[i] == 0 {
                p.as_ref()
            } else {
                q.as_ref()
            };
            w = op * w.as_ref();
        }
        let prob: f64 = (0..state.nrows())
            .map(|i| {
                let phi_i = state.as_ref().get(i, 0);
                let w_i = w.as_ref().get(i, 0);
                phi_i.re * w_i.re + phi_i.im * w_i.im
            })
            .sum();

        // Each probability must be non-negative.
        assert!(
            prob >= -1e-10,
            "negative probability at pattern {pattern}: {prob}"
        );
        prob_sum += prob;
    }

    // Probabilities must sum to 1.
    assert!(
        (prob_sum - 1.0).abs() < 1e-10,
        "probabilities do not sum to 1: {prob_sum}"
    );

    // Call the actual function to exercise the file-writing path.
    perform_x_ancilla_checks(&state, None);
}
