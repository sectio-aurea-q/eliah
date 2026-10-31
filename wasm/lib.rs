//! ELIAH for the browser (wasm32-unknown-unknown), with plain number exports and no JavaScript glue.
//!
//! `kappa` and `kappa_si` return 0 (ok), 1 (invalid input) or 2 (no convergence or overflow), like
//! the C interface, and write [value, ln_value, quad_err, tail_bound, cos_branch, evals] as six f64
//! to the buffer at `result_ptr()`. The C functions of the library (`eliah_log_transmission`,
//! `eliah_transmission_si`, ...) are exported as well.
//!
//!     cargo build --release --target wasm32-unknown-unknown
//!     node check_wasm.mjs target/wasm32-unknown-unknown/release/eliah_wasm.wasm ../conformance/eckart_certified.csv

use eliah::{EckartError, Kappa};

static mut OUT: [f64; 6] = [f64::NAN; 6];

fn store(r: Result<Kappa, EckartError>) -> i32 {
    let (code, v) = match r {
        Ok(k) => (0, [k.value, k.ln_value, k.quad_err, k.tail_bound, if k.cos_branch { 1.0 } else { 0.0 }, k.evals as f64]),
        Err(EckartError::NoConvergence) => (2, [f64::NAN; 6]),
        Err(_) => (1, [f64::NAN; 6]),
    };
    unsafe { core::ptr::addr_of_mut!(OUT).write(v) };
    code
}

/// Address of the six-number result buffer in linear memory.
#[no_mangle]
pub extern "C" fn result_ptr() -> *const f64 {
    core::ptr::addr_of!(OUT) as *const f64
}

/// κ in reduced units: barriers v1, v2 in k_B T, u = hν/(k_B T).
#[no_mangle]
pub extern "C" fn kappa(v1: f64, v2: f64, u: f64) -> i32 {
    store(eliah::kappa(v1, v2, u))
}

/// κ in SI units: barriers in J/mol, |imaginary frequency| in cm⁻¹, T in K.
#[no_mangle]
pub extern "C" fn kappa_si(v1_j_mol: f64, v2_j_mol: f64, nu_cm: f64, t_k: f64) -> i32 {
    store(eliah::kappa_si(v1_j_mol, v2_j_mol, nu_cm, t_k))
}
