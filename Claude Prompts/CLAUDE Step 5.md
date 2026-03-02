## Recap of things until now
So far, you have helped me write Rust code to simulate the surface code with coherent errors and save results in `.csv` file. Our final step is to replicate this in Python.

## Instructions
1. Make a`test_claude/python/` folder
2. Within the `test_claude/python/` folder, write a Python version of all the rust code in `test_claude/src/`.
3. For every function in rust, write an equivalent function in python with the same name. Instead of `faer-rs` in rust, use `numpy` for Python matrices. For running the Python code, **Always use the `arch` micromamba environment**. `Numpy` will already be available in the `arch` environment.
4. **DO NOT install any new python packages.** All python packages you should need should already be available in the `arch` environment. If this is not the case, alert me and check with me on what needs to be installed.
5. For every unit test in rust, write an equivalent function in python with the same name using `pytest`. If you use the `arch` environment as instructed, `pytest` will already be available.
6. Finally, write a `main.py` similar to to `main.rs` and place in within the `test_claude/python/` folder. Check that both Rust and Python versions of this simulationsgive the same values in their output `.csv` files.