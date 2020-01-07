use crate::util::sparse::Sparse;

mod gauss_siedel;
mod conjugate_gradient;

pub use self::gauss_siedel::*;
pub use self::conjugate_gradient::*;

pub enum LinearSolver {
    GaussSiedel {
        iterations: usize
    },
    ConjugateGradient {
        auxiliary: Vec<f64>,
        search: Vec<f64>,
        preconditioner: Vec<f64>,
        a: Sparse,
        iterations: usize
    }
}

impl LinearSolver {
    pub fn solve(&mut self, pressure: &mut Vec<f64>, residual: &mut Vec<f64>, cell: &Vec<u8>, fluid_density: f64, timestep: f64, cell_size: f64, rows: usize, columns: usize) {
        match self {
            LinearSolver::GaussSiedel {
                iterations
            } => {
                gauss_siedel(pressure, residual, fluid_density, timestep, cell_size, rows, columns, *iterations)
            }
            LinearSolver::ConjugateGradient {
                auxiliary,
                search,
                preconditioner,
                a,
                iterations
            } => {
                conjugate_gradient(pressure, residual, auxiliary, search, preconditioner, a, cell, fluid_density, timestep, cell_size, rows, columns, *iterations)
            }
        }
    }
}