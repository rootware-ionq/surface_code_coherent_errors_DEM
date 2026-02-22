## Picking up from Step 2
In step 2, we prepared the `logical_plus_state` (denoted mathematically as $\ket{+}_L$). Our goals now are :
- normalizing this state
- applying  coherent errors on qubits
- applying another round of syndrome measurement
## First goal: Normalizing the state
When we created the`logical_plus_state`  $\ket{+}_L$ state, we did not check if it was properly normalized. Correct normalization is essential for conserving probabilities and ensuring the sum to 1. We are going to define a function `fn renormalize_state(&mut Mat<c64>)-> Mat<c64>`. This function will take in a given $\ket{+}_L$ state, calculate its norm, and renormlize $\ket{+}_L$ . Mathematically,  this means that this function takes in $\ket{+}_L$ and
1. For each entry $c_n$ in $\ket{+}_L$ , calculate $|c_n|^2$. 
2. Calculate the norm $\mathcal{N}$  by summing $\mathcal{N} =\sum_n |c_n|^2$
3. Renormalize $\ket{+}_L$ by $\ket{+}_L \leftarrow  \frac{\ket{+}_L}{\sqrt{\mathcal{N}}}$ 
4. Return the normalized $\ket{+}_L$

To define a unit test for `renormalize_state(&mut Mat<c64>)`, give it a test state with 3 complex entries $\ket{\psi} = \frac{(1+i, i, 1)}{\sqrt{1}}$ and renormalize it. If this test is correct, it will note that $\mathcal{N} =6$ and the renormalized state returned by the function should be $\ket{\psi} = \frac{(1+i, i, 1)}{\sqrt{6}}$.

Finally, in `main`, we use `renormalize_state` to renormalize the `logical_plus_state`. We'll use this renormalized state in the next step.
## Applying Coherent Errors
We now are going to apply coherent errors to all the qubits. To do this, we will first define `get_coherent_error_matrix(theta : f64, n_qubits : usize) -> Mat<c64>`. This function will do the following:
1. First, use the given angle $\theta$ to define a single qubit coherent error matrix using the Pauli identity $I$ and Pauli $Z$ matrix
		$$ U = \cos \theta I - i Z \sin \theta$$
2. Second, we now construct $\mathcal{U}$, the coherent error matrix for `n_qubits` , by Kronecker producting $U$ `n_qubit` times:
		$$\mathcal{U} = \otimes^{n} U$$
		where $\otimes$ denotes the Kronecker product operation
3. Return the complex matrix $\mathcal{U}$

To define a unit test for function `get_coherent_error_matrix`, we will need to check:
- Is the matrix $\mathcal{U}$ Hermitian? You can test this by using $\theta = 0.1\pi$, `n_qubits = 3`.
- Is the matrix $\mathcal{U}$ the right dimensions? It should be an $2^{N} \times 2^N$ matrix, where $N=$`n_qubits`.
- For $\theta = 0$ , is $\mathcal{U}$ equal to the identity matrix?
- For $\theta = -\pi/2$, is $\mathcal{U}$ the same as Pauli $Z$ operating on all qubits?

We now define a new function that takes in a state, and applies the coherent error on all qubits. This function has the signature `fn apply_coherent_error(psi: Mat<c64>, theta: f64) -> Mat<c64>`. To do this, it
- looks at the size of `psi` to see how many qubits it corresponds to. The length of `psi` should be $2^N$ , where $N$ is the number of qubits
- calculates the coherent error matrix $\mathcal{U}$ using `get_coherent_error_matrix` for $N$ qubits using the given $\theta=$`theta`
- Multiplies $\mathcal{U}$ to `psi`= $\ket{\psi}$ , and returns the new matrix. 
Write tests for this function.

Finally, in `main`, we take the `logical_plus_state`, define $\theta = 0.1\pi$ , and apply the coherent error $\mathcal{U}$ on the `logical_plus_state`. We call this the `coherent_error_state`

## In main, printing out
In `main`, for every $i$-th entry in `coherent_error_state`, print on separate lines:
	- $i$ as a 9 bit binary number `\t`  complex norm square of `coherent_error_state[i]` 