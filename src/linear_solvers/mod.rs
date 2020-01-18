use crate::util::sparse::Sparse;

mod gauss_siedel;
mod conjugate_gradient;

pub use self::gauss_siedel::*;
pub use self::conjugate_gradient::*;
use crate::util::fluid_quantity::FluidQuantity;

pub enum LinearSolver {
    GaussSiedel,
    ConjugateGradient
}

impl LinearSolver {
    pub fn solve(&mut self,
                 pressure: &mut Vec<f64>,
                 residual: &mut Vec<f64>,
                 auxiliary: &mut Vec<f64>,
                 search: &mut Vec<f64>,
                 preconditioner: &mut Vec<f64>,
                 a: &mut Sparse,
                 cell: &Vec<u8>,
                 fluid_density: f64,
                 timestep: f64,
                 cell_size: f64,
                 rows: usize,
                 columns: usize,
                 iterations: usize,
                 u_velocity: &FluidQuantity,
                 v_velocity: &FluidQuantity) {
        match self {
            LinearSolver::GaussSiedel => {
                gauss_siedel(pressure, residual, fluid_density, timestep, cell_size, rows, columns, iterations)
            }
            LinearSolver::ConjugateGradient => {
                conjugate_gradient(pressure, residual, auxiliary, search, preconditioner, a, cell, fluid_density, timestep, cell_size, rows, columns, iterations, u_velocity, v_velocity)
            }
        }
    }
}