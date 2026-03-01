use surface_code_dem::errors::{apply_coherent_error, apply_coherent_errors_from_list};
use surface_code_dem::measurements::perform_x_ancilla_checks;
use surface_code_dem::state::{
    get_logical_plus_state, get_modified_logical_state, renormalize_state,
};

fn main() {
    println!("First test: uniform coherent errors adding constructively");
    let mut logical = get_logical_plus_state();
    let logical_plus_state = renormalize_state(&mut logical);

    let theta = 0.1 * std::f64::consts::PI;
    let coherent_error_state = apply_coherent_error(logical_plus_state, theta);

    perform_x_ancilla_checks(&coherent_error_state, None);

    println!("Second test: uniform coherent errors but some add destructively.");

    let mut logical = get_modified_logical_state();
    let logical_plus_state = renormalize_state(&mut logical);

    let theta = 0.1 * std::f64::consts::PI;
    let coherent_error_state = apply_coherent_error(logical_plus_state, theta);

    perform_x_ancilla_checks(&coherent_error_state, Some("x_ancilla_probs_2.csv"));

    println!("Third test: Non-uniform coherent errors, some additive, some destructive.");

    let mut logical = get_logical_plus_state();
    let logical_plus_state = renormalize_state(&mut logical);
    let theta_list: Vec<f64> = (0..9)
        .map(|x| {
            if x != 5 {
                0.1 * std::f64::consts::PI
            } else {
                -0.1 * std::f64::consts::PI
            }
        })
        .collect();
    let coherent_error_state = apply_coherent_errors_from_list(logical_plus_state, theta_list);

    perform_x_ancilla_checks(&coherent_error_state, Some("x_ancilla_probs_3.csv"));
}
