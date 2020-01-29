use crate::linear_solvers::LinearSolver;
use crate::integration::Integration;
use crate::advection::Advection;
use crate::util::fluid_quantity::FluidQuantity;
use crate::interpolation::Interpolation;
use crate::boundary::SolidBody;
use crate::util::helper::{max, min};
use crate::util::sparse::Sparse;

use std::time::Instant;

pub struct FluidSolver {
    pub u_velocity:     FluidQuantity,
    pub v_velocity:     FluidQuantity,
    pub density:        FluidQuantity,
    pub rows:           usize,
    pub columns:        usize,
    pressure:           Vec<f64>,
    residual:           Vec<f64>,
    auxiliary:          Vec<f64>,
    search:             Vec<f64>,
    preconditioner:     Vec<f64>,
    a:                  Sparse,
    iterations:         usize,
    timestep:           f64,
    cell_size:          f64,
    fluid_density:      f64,
    linear_solver:      LinearSolver,
    integration:        Integration,
    interpolation:      Interpolation,
    advection:          Advection,
    bodies:             Vec<SolidBody>
}

impl FluidSolver {
    // Creates a new FluidSolver. x_velocity has one more column, while y_velocity has 1 more row
    pub fn new(rows: usize, columns: usize, timestep: f64, cell_size: f64, fluid_density: f64, bodies: Vec<SolidBody>) -> FluidSolver {
        FluidSolver {
            u_velocity:     FluidQuantity::new(rows, columns + 1, 0.0, 0.5, cell_size),
            v_velocity:     FluidQuantity::new(rows + 1, columns, 0.5, 0.0, cell_size),
            density:        FluidQuantity::new(rows, columns, 0.5, 0.5, cell_size),
            rows,
            columns,
            pressure:       vec![0.0; rows * columns],
            residual:       vec![0.0; rows * columns],
            auxiliary:      vec![0.0; rows * columns],
            search:         vec![0.0; rows * columns],
            preconditioner: vec![0.0; rows * columns],
            a:              Sparse::new(rows * columns),
            iterations:     600,
            timestep,
            cell_size,
            fluid_density,
            linear_solver:  LinearSolver::GaussSiedel,
            integration:    Integration::BogackiShampine,
            interpolation:  Interpolation::BiLinear,
            advection:      Advection::SemiLagrangian,
            bodies
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
            for column in 0..self.columns {
                if self.density.cell_at(row, column) == 1 {
                    let body = &self.bodies[self.density.body_at(row, column) as usize];

                    *self.u_velocity.at_mut(row, column) = body.velocity_x(column as f64 * self.cell_size, (row as f64 + 0.5) * self.cell_size);
                    *self.v_velocity.at_mut(row, column) = body.velocity_y((column as f64 + 0.5) * self.cell_size, row as f64 * self.cell_size);
                    *self.u_velocity.at_mut(row, column + 1) = body.velocity_x((column as f64 + 1.0) * self.cell_size, (row as f64 + 0.5) * self.cell_size);
                    *self.v_velocity.at_mut(row + 1, column) = body.velocity_y((column as f64 + 0.5) * self.cell_size, (row as f64 + 1.0) * self.cell_size);
                }
            }
        }

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
                if self.density.cell_at(row, column) == 0 {
                    // Get x and x+1 velocities
                    let u1 = self.u_velocity.at(row, column) * self.u_velocity.volume_at(row, column);
                    let u2 = self.u_velocity.at(row, column + 1) * self.u_velocity.volume_at(row, column + 1);

                    // Get y and y+1 velocities
                    let v1 = self.v_velocity.at(row, column) * self.v_velocity.volume_at(row, column);
                    let v2 = self.v_velocity.at(row + 1, column) * self.v_velocity.volume_at(row + 1, column);

                    // Factor in cell scale and invert for solving
                    self.residual[row * self.columns + column] = -1.0 * (u2 - u1 + v2 - v1) / self.cell_size;

                    if !self.bodies.is_empty() {
                        let vol = self.density.volume_at(row, column);
                        let index = row * self.columns + column;

                        if column > 0 {
                            self.residual[index] -= (self.u_velocity.volume_at(row, column) - vol) * self.bodies[self.density.body_at(row, column - 1) as usize].velocity_x(column as f64 * self.cell_size, (row as f64 + 0.5) * self.cell_size);
                        }

                        if row > 0 {
                            self.residual[index] -= (self.v_velocity.volume_at(row, column) - vol) * self.bodies[self.density.body_at(row - 1, column) as usize].velocity_y((column as f64 + 0.5) * self.cell_size, row as f64 * self.cell_size);
                        }

                        if column < self.columns - 1 {
                            self.residual[index] += (self.u_velocity.volume_at(row, column + 1) - vol) * self.bodies[self.density.body_at(row, column + 1) as usize].velocity_x((column as f64 + 1.0) * self.cell_size, (row as f64 + 0.5) * self.cell_size);
                        }

                        if row < self.rows - 1 {
                            self.residual[index] += (self.v_velocity.volume_at(row + 1, column) - vol) * self.bodies[self.density.body_at(row + 1, column) as usize].velocity_y((column as f64 + 0.5) * self.cell_size, (row as f64 + 1.0) * self.cell_size);
                        }
                    }
                } else {
                    self.residual[row * self.columns + column] = 0.0;
                }

            }
        }
    }

    // Solves pressure array based on divergence, passed to linear solver function
    fn solve_pressure(&mut self) {
        let pressure_time = Instant::now();
        self.linear_solver.solve(&mut self.pressure,
                                 &mut self.residual,
                                 &mut self.auxiliary,
                                 &mut self.search,
                                 &mut self.preconditioner,
                                 &mut self.a,
                                 &self.density.cell,
                                 self.fluid_density,
                                 self.timestep,
                                 self.cell_size,
                                 self.rows,
                                 self.columns,
                                 self.iterations,
                                 &self.u_velocity,
                                 &self.v_velocity);
        print!("Linear Solve: {} ms, ", pressure_time.elapsed().as_millis())
    }

    // Applies computed pressure field to the xy velocity vector field
    fn apply_pressure(&mut self) {
        let scale = self.timestep / (self.fluid_density * self.cell_size);

        for row in 0..self.rows {
            for column in 0..self.columns {
                if self.density.cell_at(row, column) == 0 {
                    let element = row * self.columns + column;
                    *self.u_velocity.at_mut(row, column) -= scale * self.pressure[element];
                    *self.u_velocity.at_mut(row, column + 1) += scale * self.pressure[element];
                    *self.v_velocity.at_mut(row, column) -= scale * self.pressure[element];
                    *self.v_velocity.at_mut(row + 1, column) += scale * self.pressure[element];
                }
            }
        }
    }

    // Projection method implements each step of the calculation
    fn project(&mut self) {
        self.calculate_residual();
        self.solve_pressure();
        self.apply_pressure();

        self.u_velocity.extrapolate();
        self.v_velocity.extrapolate();
        self.density.extrapolate();
    }

    // Advection method moves density scalar field through velocity vector field to produce output
    fn advect(&mut self) {
        let advect_time = Instant::now();
        self.advection.advect(&mut self.u_velocity, &mut self.v_velocity, &mut self.density, self.timestep, &self.interpolation, &self.integration);
        print!("Advection: {} ms, ", advect_time.elapsed().as_millis())
    }

    // Produces the next frame of the simulation by projecting then advecting
    pub fn update(&mut self) {
        let total_time = Instant::now();
        for body in &mut self.bodies {
            body.update(self.timestep);
        }

        self.u_velocity.fill_solid_fields(&self.bodies);
        self.v_velocity.fill_solid_fields(&self.bodies);
        self.density.fill_solid_fields(&self.bodies);

        self.set_boundaries();
        self.project();
        self.set_boundaries();
        self.advect();
        println!("Total: {} ms", total_time.elapsed().as_millis())
    }

    pub fn add_inflow(&mut self, x: f64, y: f64, width: f64, height: f64, density: f64, u_velocity: f64, v_velocity: f64) {
        self.density.add_inflow(x, y, x + width, y + height, density);
        self.u_velocity.add_inflow(x, y, x + width, y + height, u_velocity);
        self.v_velocity.add_inflow(x, y, x + width, y + height, v_velocity);
    }

    // Basic function to convert density_src array into an image buffer
    pub fn to_image(&self, max_density: f64, buffer: &mut Vec<u8>) {
        for i in 0..(self.rows * self.columns) {
            let mut shade = ((max_density - self.density.src[i]) * 255.0 / max_density) as u8;
            shade = max(min(shade, 255), 0);

            if self.density.cell[i] == 1 {
                shade = 0;
            }

            buffer[i * 3 + 0] = shade;
            buffer[i * 3 + 1] = shade;
            buffer[i * 3 + 2] = shade;
        }
    }
}