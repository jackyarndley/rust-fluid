use crate::grid::{VectorGrid2, ScalarGrid2};
use crate::grid::vector2::FaceCenteredVectorGrid2;
use std::mem::swap;
use crate::linear_solver::LinearSolver;
use crate::integration::Integration;
use crate::interpolation::Interpolation;
use crate::advection::Advection;

pub struct Solver2 {
    pub velocity_src:     FaceCenteredVectorGrid2,
    pub velocity_dst:     FaceCenteredVectorGrid2,
    pub density_src:      ScalarGrid2,
    pub pressure:         ScalarGrid2,
    pub resolution:       usize,
    pub fluid_density:    f64,
    rows:                 usize,
    columns:                 usize,
    linear_solver:        LinearSolver,
    integration:          Integration,
    interpolation:        Interpolation,
    advection:            Advection
}

impl Solver2 {
    pub fn new(rows: usize, columns: usize, cell_size: f64, fluid_density: f64, linear_solver: LinearSolver, integration: Intergration, interpolation: Interpolation, advection: Advection) -> Solver2 {
        Solver2 {
            velocity_src: FaceCenteredVectorGrid2::new(rows, columns),
            velocity_dst: FaceCenteredVectorGrid2::new(rows, columns),
            density_src:  ScalarGrid2::new(rows, columns),
            pressure:     ScalarGrid2::new(rows, columns),
            resolution:   0,
            fluid_density,
            rows,
            columns,
            linear_solver,
            integration,
            interpolation,
            advection,
        }
    }

    fn set_boundaries(&mut self) {





    }

    fn project(&mut self) {
        self.calculate_residual();
        self.solve_pressure();
        self.apply_pressure();
    }

    fn swap(&mut self) {
        swap(&self.velocity_src, &self.velocity_dst)
    }

    fn advect(&mut self) {
        self.
    }

    pub fn solve(&mut self, time_step: f64) {
        self.set_boundaries();
        self.project();
        self.set_boundaries();
        self.advect();
    }

    pub fn generate_image(&self, buffer: &mut Vec<u8>) {
        for i in 0..self.density_src.len {
            let shade = (self.density_src * 255.0 / 3.0) as u8;

            buffer[i * 4 + 0] = shade;
            buffer[i * 4 + 1] = shade;
            buffer[i * 4 + 2] = shade;
            buffer[i * 4 + 3] = 0xFF;
        }
    }
}