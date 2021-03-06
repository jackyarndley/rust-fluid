<p align="center">
  <img src="examples/image.png" width="100%">
</p>

# rust-fluid
Implementation of a 2-dimensional MAC (marker and cell) fluid solver in Rust. A variety of core algorithms are implemented and can be easily changed and swapped out to customise the solving method.

## Usage
Currently, the primary method of modifying the scene and parameters of the simulation is by modifying the source code. This is done using the main.rs file. An example scene is already provided.

```cargo run --release```

The output files are placed in the ```/output``` directory. This directory must be created prior to running the solver.

## Resources
- [incremental-fluids](https://github.com/tunabrain/incremental-fluids) provided comprehensive documentation and implementations of the algorithms used in this project.
