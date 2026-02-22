use surface_code_dem::errors::apply_coherent_error;
use surface_code_dem::measurements::perform_x_ancilla_checks;
use surface_code_dem::state::{get_logical_plus_state, renormalize_state};

fn main() {
    let mut logical = get_logical_plus_state();
    let logical_plus_state = renormalize_state(&mut logical);

    let theta = 0.1 * std::f64::consts::PI;
    let coherent_error_state = apply_coherent_error(logical_plus_state, theta);

    perform_x_ancilla_checks(&coherent_error_state);
}
