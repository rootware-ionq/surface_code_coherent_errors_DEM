## Projecting out the logical state
We have previously defined functionality for building the $X$ and $Z$ stabilizers in the previous step. Confirm this.

### projection operations
We now write a function to take in a stabilizer $S$(can be $X$ or $Z$ stabilizer matrix) and return the projection operator $P = \frac{ \mathbb{I} + S}{2}$ . Here, $\mathbb{I}$ is the Pauli identity $I$ on all 9 qubits. Write this function and test it. Note that $P$ should have the same size and dimensions as $S$.

We now generate projection operators for all the possible $X$ and $Z$ stabilizers, and then multiply them sequentially onto the initial state. We then record this new state as our `logical_plus_state` . Mathematically, we are making this new state $\ket{+}_L$ by

$\ket{+}_L = \left(\prod_{S\in stabilizers} \frac{\mathbb{I}+S}{2} \right)\ket{\psi_0}$

where $\ket{\psi_0}$ is the state returned by the `get_initial_state()` function for `n_qubits=9`

Print out this `logical_plus_state`
