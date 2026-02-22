## Mandatory Instructions

1. For each function you define, and check with me before editing file.
2. For each function you define, write a test, and check with me before editing file.
3. When in doubt, ask me for clarification.
4. Go as step by step as possible.
## Project Summary

**Goal**: You are about to implement a specific simulation algorithm for simulating the rotated surface code, a quantum error correcting code. You will simulate the effect and failure rate of the rotated surface code in the presence of coherent errors

## Libraries to use
1. Language: Rust
2. For matrix definitions and linear algebra: `faer-rs`
3. By default, use complex numbers as data type (`c64` in `faer-rs`).
## Step 1 : Defining building blocks

We first define functions relevant to state initialization.
1. Define a function for a single qubit being in the $\ket{+} = \frac{1}{2} (\ket{0}+ \ket{1})$ state. This function should have the definition that looks something like
		`fn plus_ket() -> Mat<c64>`
2. Now define a function named `get_initial_state()` that takes in number of qubits and defines their initial wavefunction as all of the qubits being in the $\ket{+}$ state. You can use the `kron` function in `faer-rs` for this. The `get_initial_state` function has a definition that looks like
		`fn get_initial_state( n_qubits: usize) -> Mat<c64>`

We now work towards defining the stabilizer operators for the $d=3$ rotated surface code. The $d=3$ rotated surface code has 9 data qubits only.  We label the qubits from 0 to 8, then
1. First, we define the Pauli matrices $X$,$Y$,$Z$ for 1 qubit. Each of these should be `Mat<c64>`
2. An $X$ stabilizer is defined as a product of the Pauli $X$'s operating on a set of qubits, and an identity on all the remaining ones. For example, a $X$ stabilizer could be $X_{stab} = I_0 X_1X_2I_3X_4X_5 I_6I_7I_8$. To simplify, we can write this stabilizer by recording only the qubits where we operate an $X$ , e.g. $X_{stab} \to (1,2,4,5)$ .
3. A $Z$ stabilizer is defined as a product of the Pauli $Z$'s operating on a set of qubits, and an identity on all the remaining ones. For example, a $Z$ stabilizer could be $Z_{stab} = I_0 I_1Z_2I_3I_4Z_5 I_6I_7I_8$. To simplify, we can write this stabilizer by recording only the qubits where we operate an $Z$ , e.g. $Z_{stab} \to (2,5)$ .
4. We then define `Vec<Vec<usize>>` vectors that contain the indices for  the qubits participating in the $X$ and $Z$ stabilizers in the $d=3$ rotated surface code.
	1. The $X$ stabilizers  operate on sets from the following list of qubits : $((0,1), (1,2,4,5), (3,4,6,7), (7,8))$. We can write this as a `Vec<Vec<usize>>`
	2. The $Z$ stabilizers operate on sets from the following list of qubits: $((0,1,3,4), (2,5), (3,6), (4,5,7,8))$. We can write this as a `Vec<Vec<usize>>`
5. Construct the matrix representation of the $X$ stabilizers and the $Z$ stabilizers by Kroneck producting the identity $I$ ,$X$, and $Z$ for each stabilizer, using the indices we recorded in `Vec<Vec<usize>>`