use crate::linear_solvers::LinearSolver;
use crate::integration::Integration;
use crate::advection::Advection;
use crate::util::fluid_quantity::FluidQuantity;
use crate::interpolation::Interpolation;

pub struct FluidSolver {
    pub u_velocity: FluidQuantity,
    pub v_velocity: FluidQuantity,
    pub density: FluidQuantity,
    pub rows:           usize,
    pub columns:        usize,
    pressure:           Vec<f64>,
    residual:           Vec<f64>,
    timestep:           f64,
    cell_size:          f64,
    fluid_density:      f64,
    linear_solver:      LinearSolver,
    integration:        Integration,
    interpolation:      Interpolation,
    advection:          Advection,
}

impl FluidSolver {
    // Creates a new FluidSolver. x_velocity has one more column, while y_velocity has 1 more row
    pub fn new(rows: usize, columns: usize, timestep: f64, cell_size: f64, fluid_density: f64) -> FluidSolver {
        FluidSolver {
            u_velocity:     FluidQuantity::new(rows, columns + 1, 0.0, 0.5, cell_size),
            v_velocity:     FluidQuantity::new(rows + 1, columns, 0.5, 0.0, cell_size),
            density:        FluidQuantity::new(rows, columns, 0.5, 0.5, cell_size),
            rows,
            columns,
            pressure:       vec![0.0; rows * columns],
            residual:       vec![0.0; rows * columns],
            timestep,
            cell_size,
            fluid_density,
            linear_solver:  LinearSolver::GaussSiedel {
                iterations: 600
            },
            integration:    Integration::BogackiShampine,
            interpolation:  Interpolation::BiLinear,
            advection:      Advection::SemiLagrangian,
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

        // TODO: Set all bordering solids to the solid velocity

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
                    let u1 = self.u_velocity.at(row, column);
                    let u2 = self.u_velocity.at(row, column + 1);

                    // Get y and y+1 velocities
                    let v1 = self.v_velocity.at(row, column);
                    let v2 = self.v_velocity.at(row + 1, column);

                    // Factor in cell scale and invert for solving
                    self.residual[row * self.columns + column] = -1.0 * (u2 - u1 + v2 - v1) / self.cell_size;
                } else {
                    self.residual[row * self.columns + column] = 0.0;
                }

            }
        }
    }

    // Solves pressure array based on divergence, passed to linear solver function
    fn solve_pressure(&mut self) {
        self.linear_solver.solve(&mut self.pressure, &mut self.residual, self.fluid_density, self.timestep, self.cell_size, self.rows, self.columns);
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
    }

    // Advection method moves density scalar field through velocity vector field to produce output
    fn advect(&mut self) {
        self.advection.advect(&mut self.u_velocity, &mut self.v_velocity, &mut self.density, self.timestep, &self.interpolation, &self.integration);
    }

    // Produces the next frame of the simulation by projecting then advecting
    pub fn update(&mut self) {
        self.set_boundaries();
        self.project();
        self.set_boundaries();
        self.advect();
    }

    // TODO: Add an inflow function

    // Basic function to convert density_src array into an image buffer
    pub fn to_image(&self, max_density: f64, buffer: &mut Vec<u8>) {
        for i in 0..(self.rows * self.columns) {
            let shade = (self.density.src[i] * 255.99 / max_density).trunc() as u8;

            buffer[i * 4 + 0] = shade;
            buffer[i * 4 + 1] = shade;
            buffer[i * 4 + 2] = shade;
            buffer[i * 4 + 3] = 0xFF; // Alpha is 255
        }
    }
}