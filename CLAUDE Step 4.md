## Picking up from Step 3
In Step 3, we prepared a logical state $\ket{+}_L$ , normalized it, and then applied coherent errors to it.

## Making stabilizer/ancilla measurements
### Writing "anti-projection" operator
Back in Step 2, we defined a function `pub fn projection_operator(s: &Mat<c64>) -> Mat<c64>` that takes in a stabilizer $S$ and returns the projection $P = \frac{ \mathbb{I} + S}{2}$. $P$ projects a state into the +1 eigenstate of $S$.

I added a `pub fn anti_projection_operator(s: &Mat<c64>) -> Mat<c64>` that takes in a stabilizer $S$ and returns the antiprojection $Q = \frac{ \mathbb{I} - S}{2}$. $Q$ projects a state into the -1 eigenstate of $S$.

Write tests for this function similar to how you wrote tests for projection operators

### Calculating probabilities for different $X$ stabilizer measurements
Our goal now is to measure the probabilities for measuring a given set of $X$ ancilla measurements on the state `coherent_error_state`. We're going to define a function `perform_x_ancilla_checks(state: Mat<c64>)` that does this. This function will take in a given state, and do the following.

- Recall first that we have `build_x_stabilizers()` function that constructs a vector of $X$ stabilizers. We call this `vec_x_stab` (mathematically, $\vec{X}_{stab}$)
- Mathematically, we define ancilla measurements as $A=(a_0,a_1,a_2,a_3)$ . $a_i$ can be $0$ or $1$. If $a_i$ is 0, it means measuring the +1 eigenvalue for the $i-$th $X$ stabilizer in $\vec{X}_{stab}$ . If $a_i=1$ ,  it means measuring the -1 eigenvalue for the $i-$th $X$ stabilizer in $\vec{X}_{stab}$ .
- Iterate over all possible ancilla measurements $A=(a_0,a_1,a_2,a_3)$ , starting from $A=(0,0,0,0)$ to $(1,1,1,1)$ and:
	- we're going to construct the ancilla measurement operator $\mathcal{O}=O_0O_1O_2O_3$  for the current ancilla.
	- If $a_i=0$, $O_i=P_i$ where $P_i$ is the projection operation for the $i-th$ $X$ stabilizer in $\vec{X}_{stab}$ 
	- If $a_i=1$, $O_i=Q_i$ where $Q_i$ is the anti-projection operation for the $i-th$ $X$ stabilizer in $\vec{X}_{stab}$ 
	- Denoting $\ket{\phi}$ as the `coherent_error_state` state, calculate the probability for the current ancilla $A$ as $p_A =\braket{\phi|\mathcal{O}|\phi}$ .
	- Record the the tuple $A$ and corresponding probability $p_A$ 
	- Move to the next ancilla
- print to a file all ancilla tuples $A$ and their corresponding probabilities $p_A$ , in a format that will be easier to read in Python NumPy code scripts in future.

Finally, in `main` , add code that calls `perform_x_ancilla_checks` on `coherent_error_state`.
### Writing Tests
Write tests for the above functions, and add documentation.