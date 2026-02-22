# $d=3 \quad X-$ memory rotated surface code 
This repo implements Rust code to simulate the experiments for the $d=3$ rotated surface code with $X-$ memory checks from this [paper by Evangelia Takou et al.](https://arxiv.org/abs/2510.23797). The algorithm in this code is shared by the paper's author through private correspondence with Shah Saad Alam. I (Shah) implemented it both as an exercise to better understand the paper for a journal club, understand coherent errors in DEMs as well as learn how to use Claude better 🤖

For my journal club talk on this paper, see recordings under Architecture Journal Club schedule.

**Acknowledgements**: Massive thanks to the author E. Takou for explaining their paper to me. This code is more of a learning exercise I did based on conversations with the author.

### Building the repo and simulating circuit with coherent errors

- Install Rust and Cargo using `rustup`
- Within the directory, run `cargo run --release`

This should generate the `x_ancilla_probs.csv` file.

### Building the DEM
Use the notebook `building the DEM.ipynb` to read the `x_ancilla_probs.csv` file, and build the DEM by inferring the probabilities $p_{ij}$. The notebook also translates these probabilities into a effective coherent angle $\tilde{theta}_{ij}$ using $p_{ij}=\sin^2 \tilde{\theta}_{ij}$. Note that for the example `main.rs` and `.csv` provided in this repo, $\theta=0.1\pi$ for the coherent errors, so we expect most $\theta_{ij}$ to be equal to $\theta$, except for the weight 4 checks where 2 data qubits contribute to one DEM edge. For those, we expect $\theta_{ij} = 2\theta$.

### Assorted Notes
For some math and derivations I did on why, see Notes.md.
