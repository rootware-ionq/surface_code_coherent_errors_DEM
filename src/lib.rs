//! Rotated surface code simulation library.
//!
//! Provides building blocks for simulating the d=3 rotated surface code
//! under coherent Z-rotation errors.
//!
//! # Modules
//! - [`gates`] — Pauli matrices and identity operators
//! - [`stabilizers`] — Stabilizer construction and projection operators
//! - [`state`] — State initialization, logical state preparation, normalization
//! - [`errors`] — Coherent error unitaries

pub mod errors;
pub mod gates;
pub mod measurements;
pub mod stabilizers;
pub mod state;
