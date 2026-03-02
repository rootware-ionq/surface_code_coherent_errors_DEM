# $d=3 \quad X-$ memory rotated surface code 
This repo implements Rust code to simulate the experiments for the $d=3$ rotated surface code with $X-$ memory checks from this [paper by Evangelia Takou et al.](https://arxiv.org/abs/2510.23797). The algorithm in this code is from the paper's author who shared it with me via private correspondence. I mostly implemented it both as an exercise to better understand the paper for a journal club, understand coherent errors in DEMs as well as learn how to use Claude better 🤖

**Acknowledgements**: Massive thanks to the author E. Takou for explaining their paper to me. This code is more of a learning exercise I did based on conversations with the author.

### `[Rust]` : Building the repo and simulating circuit with coherent errors
This is the primary code used to generate data files for the surface code. I also added some additional functionality for custom experiments that is not present in the Python code (see below). To set up:

- Install Rust and Cargo using `rustup`
- Within the directory, run `cargo run --release`
- The Rust source code is within the `/src/` folder, with tests in `/tests/`.

This should generate the `x_ancilla_probs.csv` file.

### `[Python]` : Building the repo and simulating circuit with coherent errors
The Python code for the surface code was translated from Rust to Python using Claude.

- Have a Python environment with `Numpy` and `pytest` installed.
- All Python code and tests are under `/python/` folder.

### Building the DEM

Use the notebook in `analysis notebooks/First step : building the DEM.ipynb` as an entry point. This notebook I wrote reads the `x_ancilla_probs.csv` file, and builds the DEM by inferring the probabilities $p_{ij}$. The notebook also translates these probabilities into a effective coherent angle $`\tilde{\theta}_{ij}`$ using $`p_{ij}=\sin^2 \tilde{\theta}_{ij}`$. Note that for the example `main.rs` and `.csv` provided in this repo, $\theta=0.1\pi$ for the coherent errors, so we expect most $`\theta_{ij}`$ to be equal to $\theta$, except for the weight 4 checks where 2 data qubits contribute to one DEM edge. For those, we expect $`\theta_{ij} = 2\theta`$.

### Math Notes
For some math and derivations I did on why, see `Math Notes.md`.
