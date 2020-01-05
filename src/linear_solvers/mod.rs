use crate::util::Sparse;

mod gauss_siedel;
mod conjugate_gradient;

pub use self::gauss_siedel::*;
pub use self::conjugate_gradient::*;

pub enum LinearSolver {
    Empty,
    GaussSiedel {
        limit: usize
    },
    ConjugateGradient {
        auxiliary: Vec<f64>,
        search: Vec<f64>,
        preconditioner: Vec<f64>,
        a: Sparse,
        limit: usize
    }
}

impl LinearSolver {
    pub fn run(&mut self, pressure: &mut Vec<f64>, residual: &mut Vec<f64>, fluid_density: f64, dt: f64, dx: f64, rows: usize, columns: usize) {
        match self {
            LinearSolver::Empty => {

            }
            LinearSolver::GaussSiedel {
                limit
            } => {
                gauss_siedel(pressure, residual, fluid_density, dt, dx, rows, columns, *limit)
            }
            LinearSolver::ConjugateGradient {
                auxiliary,
                search,
                preconditioner,
                a,
                limit
            } => {
                conjugate_gradient(pressure, residual, auxiliary, search, preconditioner, a, fluid_density, dt, dx, rows, columns, *limit)
            }
        }
    }
}