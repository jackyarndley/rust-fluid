# rust-fluid
## Overview
Implementation of a 2-dimensional MAC (marker and cell) fluid solver in Rust. A variety of different implementations of the code algorithms are presented and easily changed to customise the solving method.

## Usage
Currently, the primary method of modifying the scene and parameters of the simulation is by modifying the source code. This is done using the main.rs file. An example scene is already provided.

```cargo run --release```

The output files are placed in the ```/output``` directory. This directory must be created prior to running the solver.

## Structure
This project has been structured to keep all the seperate all of the different components.

## Resources
[incremental-fluids](https://github.com/tunabrain/incremental-fluids) provided comprehensive documentation and implementations of the algorithms used in this project.
