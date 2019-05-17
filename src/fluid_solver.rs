use crate::util::{LinearSolver, Integration, Interpolation, Advection};
use crate::util::{Field, Sparse};
use std::mem::swap;
use crate::linear_solvers;
use crate::integration;
use crate::interpolation;
use crate::advection;

pub struct FluidSolver {
    pub u_velocity:     Field,
    pub u_velocity_dst: Field,
    pub v_velocity:     Field,
    pub v_velocity_dst: Field,
    pub density:        Field,
    pub density_dst:    Field,
    pub rows:           usize,
    pub columns:        usize,
    pressure:           Vec<f64>,
    residual:           Vec<f64>,
    auxiliary:          Vec<f64>,
    search:             Vec<f64>,
    preconditioner:     Vec<f64>,
    a:                  Sparse,
    dt:                 f64,
    dx:                 f64,
    fluid_density:      f64,
    linear_solver:      LinearSolver,
    integration:        Integration,
    interpolation:      Interpolation,
    advection:          Advection,
}

impl FluidSolver {
    // Creates a new FluidSolver. x_velocity has one more column, while y_velocity has 1 more row
    pub fn new(rows: usize, columns: usize, dt: f64, dx: f64, fluid_density: f64) -> FluidSolver {
        FluidSolver {
            u_velocity:     Field::new(rows, columns + 1, 0.0, 0.5),
            u_velocity_dst: Field::new(rows, columns + 1, 0.0, 0.5),
            v_velocity:     Field::new(rows + 1, columns, 0.5, 0.0),
            v_velocity_dst: Field::new(rows + 1, columns, 0.5, 0.0),
            density:        Field::new(rows, columns, 0.5, 0.5),
            density_dst:    Field::new(rows, columns, 0.5, 0.5),
            rows,
            columns,
            pressure:       vec![0.0; rows * columns],
            residual:       vec![0.0; rows * columns],
            auxiliary:      vec![0.0; rows * columns],
            search:         vec![0.0; rows * columns],
            preconditioner: vec![0.0; rows * columns],
            a:              Sparse::new(rows * columns),
            dt,
            dx,
            fluid_density,
            linear_solver:  linear_solvers::empty,
            integration:    integration::empty,
            interpolation:  interpolation::empty,
            advection:      advection::empty,
        }
    }

    // Sets the linear solver used in the simulation
    pub fn linear_solver(mut self, f: LinearSolver) -> Self {
        self.linear_solver = f;
        self
    }

    // Sets the integration method used in the simulation
    pub fn integration(mut self, f: Integration) -> Self {
        self.integration = f;
        self
    }

    // Sets the interpolation method used in the simulation
    pub fn interpolation(mut self, f: Interpolation) -> Self {
        self.interpolation = f;
        self
    }

    // Sets the advection method used in the simulation
    pub fn advection(mut self, f: Advection) -> Self {
        self.advection = f;
        self
    }

    // Sets boundaries of simulation by setting xy velocities at boundaries to 0
    fn set_boundaries(&mut self) {
        for row in 0..self.rows {
            *self.u_velocity.at_mut(row, 0) = 0.0;
            *self.u_velocity.at_mut(row, self.columns) = 0.0;
        }

        for column in 0..self.columns {
            *self.v_velocity.at_mut(0, column) = 0.0;
            *self.v_velocity.at_mut(self.rows, column) = 0.0;
        }
    }

    // Calculates residual vector from uv vector field
    fn calculate_residual(&mut self) {
        for row in 0..self.rows {
            for column in 0..self.columns {
                // Get x and x+1 velocities
                let u1 = self.u_velocity.at(row, column);
                let u2 = self.u_velocity.at(row, column + 1);

                // Get y and y+1 velocities
                let v1 = self.v_velocity.at(row, column);
                let v2 = self.v_velocity.at(row + 1, column);

                // Factor in cell scale and invert for solving
                self.residual[row * self.columns + column] = -1.0 * (u2 - u1 + v2 - v1) / self.dx;
            }
        }
    }

    // Solves pressure array based on divergence, passed to linear solver function
    fn solve_pressure(&mut self) {
        (self.linear_solver)(&mut self.pressure, &mut self.residual, &mut self.auxiliary, &mut self.search, &mut self.preconditioner, &mut self.a, self.fluid_density, self.dt, self.dx, self.rows, self.columns, 600);
    }

    // Applies computed pressure field to the xy velocity vector field
    fn apply_pressure(&mut self) {
        let scale = self.dt / (self.fluid_density * self.dx);

        for row in 0..self.rows {
            for column in 0..self.columns {
                let element = row * self.columns + column;
                *self.u_velocity.at_mut(row, column) -= scale * self.pressure[element];
                *self.u_velocity.at_mut(row, column + 1) += scale * self.pressure[element];
                *self.v_velocity.at_mut(row, column) -= scale * self.pressure[element];
                *self.v_velocity.at_mut(row + 1, column) += scale * self.pressure[element];
            }
        }
    }

    // Projection method implements each step of the calculation
    fn project(&mut self) {
        self.calculate_residual();
        self.solve_pressure();
        self.apply_pressure();
    }

    // Advection method moves density scalar field through velocity vector field to produce output
    fn advect(&mut self) {
        (self.advection)(&mut self.u_velocity_dst, &self.u_velocity, &self.u_velocity, &self.v_velocity, self.dt, self.dx, &self.integration);
        (self.advection)(&mut self.v_velocity_dst, &self.v_velocity, &self.u_velocity, &self.v_velocity, self.dt, self.dx, &self.integration);
        (self.advection)(&mut self.density_dst, &self.density, &self.u_velocity, &self.v_velocity, self.dt, self.dx, &self.integration);

        swap(&mut self.u_velocity.field, &mut self.u_velocity_dst.field);
        swap(&mut self.v_velocity.field, &mut self.v_velocity_dst.field);
        swap(&mut self.density.field, &mut self.density_dst.field);
    }

    // Produces the next frame of the simulation by projecting then advecting
    pub fn solve(&mut self) {
        self.set_boundaries();
        self.project();
        self.set_boundaries();
        self.advect();
    }

    // Basic function to convert density_src array into an image buffer
    pub fn to_image(&self, buffer: &mut Vec<u8>) {
        for i in 0..(self.rows * self.columns) {
            let shade: u8 = (self.density.field[i] * 255.0 / 3.0) as u8;

            buffer[i * 4 + 0] = shade;
            buffer[i * 4 + 1] = shade;
            buffer[i * 4 + 2] = shade;
            buffer[i * 4 + 3] = 0xFF; // Alpha is 255
        }
    }
}