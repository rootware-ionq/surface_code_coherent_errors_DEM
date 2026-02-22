use surface_code_dem::errors::apply_coherent_error;
use surface_code_dem::measurements::perform_x_ancilla_checks;
use surface_code_dem::state::{get_logical_plus_state, renormalize_state};

fn main() {
    let mut logical = get_logical_plus_state();
    let logical_plus_state = renormalize_state(&mut logical);

    let theta = 0.1 * std::f64::consts::PI;
    let coherent_error_state = apply_coherent_error(logical_plus_state, theta);

    println!("coherent_error_state (index [binary]\tnorm^2):");

    let mut sum = 0.0;
    for i in 0..coherent_error_state.nrows() {
        let v = coherent_error_state.as_ref().get(i, 0);
        let norm_sq = v.re * v.re + v.im * v.im;

        if norm_sq > 1e-10 {
            println!("{i:09b}\t{norm_sq:.10}");
        }
        sum += norm_sq;
    }
    println!("1.0 - sum of coeff:\t{:e}", 1.0 - sum);

    perform_x_ancilla_checks(&coherent_error_state);
}
