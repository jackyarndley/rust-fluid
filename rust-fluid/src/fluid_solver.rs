use crate::types::*;
use crate::field::Field;
use crate::linear_solvers;
use crate::integration;
use crate::interpolation;
use crate::advection;
use std::mem::swap;

pub struct FluidSolver {
    pub x_velocity_src: Field,
    pub x_velocity_dst: Field,
    pub y_velocity_src: Field,
    pub y_velocity_dst: Field,
    pub density_src:    Field,
    pub density_dst:    Field,
    pub pressure:       Field,
    pub divergence:     Field,
    pub rows:           usize,
    pub columns:        usize,
    dt:                 f64,
    dx:                 f64,
    fluid_density:      f64,
    linear_solver:      LinearSolver,
    integration:        Integration,
    interpolation:      Interpolation,
    advection:          Advection,
}

impl FluidSolver {
    pub fn new(rows: usize, columns: usize, dt: f64, dx: f64, fluid_density: f64) -> FluidSolver {
        FluidSolver {
            x_velocity_src: Field::new(rows, columns + 1, 0.0, 0.5),
            x_velocity_dst: Field::new(rows, columns + 1, 0.0, 0.5),
            y_velocity_src: Field::new(rows + 1, columns, 0.5, 0.0),
            y_velocity_dst: Field::new(rows + 1, columns, 0.5, 0.0),
            density_src:    Field::new(rows, columns, 0.5, 0.5),
            density_dst:    Field::new(rows, columns, 0.5, 0.5),
            pressure:       Field::new(rows, columns, 0.5, 0.5),
            divergence:     Field::new(rows, columns, 0.5, 0.5),
            rows,
            columns,
            dt,
            dx,
            fluid_density,
            linear_solver:  linear_solvers::empty,
            integration:    integration::empty,
            interpolation:  interpolation::empty,
            advection:      advection::empty,
        }
    }

    pub fn linear_solver(&mut self, f: LinearSolver) {
        self.linear_solver = f
    }

    pub fn integration(&mut self, f: Integration) {
        self.integration = f
    }

    pub fn interpolation(&mut self, f: Interpolation) {
        self.interpolation = f
    }

    pub fn advection(&mut self, f: Advection) {
        self.advection = f
    }

    fn calculate_divergence(&mut self) {
        for row in 0..self.rows {
            for column in 0..self.columns {
                let x_velocity1 = self.x_velocity_src.at(row, column);
                let x_velocity2 = self.x_velocity_src.at(row, column + 1);

                let y_velocity1 = self.y_velocity_src.at(row, column);
                let y_velocity2 = self.y_velocity_src.at(row + 1, column);

                *self.divergence.at_mut(row, column) = - (x_velocity2 - x_velocity1 + y_velocity2 - y_velocity1) / self.dx;
            }
        }
    }

    fn solve_pressure(&mut self) {
        (self.linear_solver)(&mut self.pressure, &self.divergence, self.fluid_density, self.dt, self.dx, 600);
    }

    fn apply_pressure(&mut self) {
        let scale = self.dt / (self.fluid_density * self.dx);

        for row in 0..self.rows {
            for column in 0..self.columns {
                *self.x_velocity_src.at_mut(row, column) -= scale * self.pressure.at(row, column);
                *self.x_velocity_src.at_mut(row, column + 1) += scale * self.pressure.at(row, column);
                *self.y_velocity_src.at_mut(row, column) -= scale * self.pressure.at(row, column);
                *self.y_velocity_src.at_mut(row + 1, column) += scale * self.pressure.at(row, column);
            }
        }
    }

    fn set_boundaries(&mut self) {
        for row in 0..self.rows {
            *self.x_velocity_src.at_mut(row, 0) = 0.0;
            *self.x_velocity_src.at_mut(row, self.columns) = 0.0;
        }

        for column in 0..self.columns {
            *self.y_velocity_src.at_mut(0, column) = 0.0;
            *self.y_velocity_src.at_mut(self.rows, column) = 0.0;
        }
    }


    fn project(&mut self) {
        self.calculate_divergence();
        self.solve_pressure();
        self.apply_pressure();
        self.set_boundaries();
    }

    fn advect(&mut self) {
        (self.advection)(&mut self.x_velocity_dst, &self.x_velocity_src, &self.y_velocity_src, self.dt, self.dx, &self.interpolation, &self.integration);
        (self.advection)(&mut self.y_velocity_dst, &self.x_velocity_src, &self.y_velocity_src, self.dt, self.dx, &self.interpolation, &self.integration);
        (self.advection)(&mut self.density_dst, &self.x_velocity_src, &self.y_velocity_src, self.dt, self.dx, &self.interpolation, &self.integration);
        self.swap();
    }

    pub fn solve(&mut self) {
        self.project();
        self.advect();
    }

    // Swaps src and dst buffers of each required field
    fn swap(&mut self) {
        swap(&mut self.x_velocity_src.field, &mut self.x_velocity_dst.field);
        swap(&mut self.y_velocity_src.field, &mut self.y_velocity_dst.field);
        swap(&mut self.density_src.field, &mut self.density_dst.field);
    }
}