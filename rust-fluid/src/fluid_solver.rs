use crate::types::*;
use crate::field::Field;
use crate::linear_solvers;
use crate::integration;
use crate::interpolation;
use crate::advection;
use std::mem::swap;
use crate::interpolation::max;
use crate::field::Vector2D;

// FluidSolver struct containing src and dst buffers and simulation details
pub struct FluidSolver {
    pub x_velocity_src: Field,
    pub x_velocity_dst: Field,
    pub y_velocity_src: Field,
    pub y_velocity_dst: Field,
    pub density_src:    Field,
    pub density_dst:    Field,
    pub p:       Vector2D,
    pub r:     Vector2D,
    pub z:              Vector2D,
    pub s:              Vector2D,
    pub pre_con:        Vector2D,
    pub a_diag:         Vector2D,
    pub a_plus_x:       Vector2D,
    pub a_plus_y:       Vector2D,
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
    // Creates a new FluidSolver. x_velocity has one more column, while y_velocity has 1 more row
    pub fn new(rows: usize, columns: usize, dt: f64, dx: f64, fluid_density: f64) -> FluidSolver {
        FluidSolver {
            x_velocity_src: Field::new(rows, columns + 1, 0.0, 0.5),
            x_velocity_dst: Field::new(rows, columns + 1, 0.0, 0.5),
            y_velocity_src: Field::new(rows + 1, columns, 0.5, 0.0),
            y_velocity_dst: Field::new(rows + 1, columns, 0.5, 0.0),
            density_src:    Field::new(rows, columns, 0.5, 0.5),
            density_dst:    Field::new(rows, columns, 0.5, 0.5),
            p:              Vector2D::new(rows, columns),
            r:              Vector2D::new(rows, columns),
            z:              Vector2D::new(rows, columns),
            s:              Vector2D::new(rows, columns),
            pre_con:        Vector2D::new(rows, columns),
            a_diag:         Vector2D::new(rows, columns),
            a_plus_x:       Vector2D::new(rows, columns),
            a_plus_y:       Vector2D::new(rows, columns),
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

    // Calculates the divergence of the xy vector field into the divergence scalar field
    fn calculate_rhs(&mut self) {
        for row in 0..self.rows {
            for column in 0..self.columns {
                // Get x and x+1 velocities
                let x_velocity1 = self.x_velocity_src.at(row, column);
                let x_velocity2 = self.x_velocity_src.at(row, column + 1);

                // Get y and y+1 velocities
                let y_velocity1 = self.y_velocity_src.at(row, column);
                let y_velocity2 = self.y_velocity_src.at(row + 1, column);

                // For fluid simulation case divergence is negative
                *self.r.at_mut(row, column) = -1.0 * (x_velocity2 - x_velocity1 + y_velocity2 - y_velocity1) / self.dx;
            }
        }
    }

    fn build_pressure_matrix(&mut self) {
        let scale = self.dt / (self.fluid_density * self.dx * self.dx);

        self.a_diag.field = vec![0.0; self.rows * self.columns];
        self.a_plus_x.field = vec![0.0; self.rows * self.columns];
        self.a_plus_y.field = vec![0.0; self.rows * self.columns];

        for y in 0..self.rows {
            for x in 0..self.columns {
                if x < (self.columns - 1) {
                    *self.a_diag.at_mut(y, x) += scale;
                    *self.a_diag.at_mut(y, x + 1) += scale;
                    *self.a_plus_x.at_mut(y, x) = -scale;
                } else {
                    *self.a_plus_x.at_mut(y, x) = 0.0;
                }

                if y < (self.rows - 1) {
                    *self.a_diag.at_mut(y, x) += scale;
                    *self.a_diag.at_mut(y + 1, x) += scale;
                    *self.a_plus_y.at_mut(y, x) = -scale;
                } else {
                    *self.a_plus_y.at_mut(y, x) = 0.0;
                }
            }
        }


    }

    fn build_preconditioner(&mut self) {
        let tau = 0.97;

        for y in 0..self.rows {
            for x in 0..self.columns {
                let mut e = self.a_diag.at(y, x);

                if x > 0 {
                    let px = self.a_plus_x.at(y, x - 1) * self.pre_con.at(y, x - 1);
                    let py = self.a_plus_y.at(y, x - 1) * self.pre_con.at(y, x - 1);
                    e -= px * px + tau * px * py;
                }

                if y > 0 {
                    let px = self.a_plus_x.at(y - 1, x) * self.pre_con.at(y - 1, x);
                    let py = self.a_plus_y.at(y - 1, x) * self.pre_con.at(y - 1, x);
                    e -= py * py + tau * px * py;
                }



                *self.pre_con.at_mut(y, x) = 1.0 / ((e + 1e-30).sqrt());
            }
        }
    }

    fn apply_preconditioner(dst: &mut Vector2D, a: &Vector2D, a_plus_x: &Vector2D, a_plus_y: &Vector2D, pre_con: &Vector2D) {
        for y in 0..a.rows {
            for x in 0..a.columns {
                let mut t = a.at(y, x);

                if x > 0 {
                    t -= a_plus_x.at(y, x - 1) * pre_con.at(y, x - 1) * dst.at(y, x - 1);
                }

                if y > 0 {
                    t -= a_plus_y.at(y - 1, x) * pre_con.at(y - 1, x) * dst.at(y - 1, x);
                }

                *dst.at_mut(y, x) = t * pre_con.at(y, x);
            }
        }

        for y in (0..a.rows).rev() {
            for x in (0..a.columns).rev() {
                let mut t = dst.at(y, x);

                if x < (a.columns - 1) {
                    t -= a_plus_x.at(y, x) * pre_con.at(y, x) * dst.at(y, x + 1);
                }

                if y < (a.rows - 1) {
                    t -= a_plus_y.at(y, x) * pre_con.at(y, x) * dst.at(y + 1, x);
                }

                *dst.at_mut(y, x) = t * pre_con.at(y, x);
            }
        }
    }

    fn dot_product(a: &Vector2D, b: &Vector2D) -> f64 {
        let mut result = 0.0;
        for y in 0..a.rows {
            for x in 0..a.columns {
                result += a.at(y, x) * b.at(y, x);
            }
        }
        result
    }

    // Multiplies pressure matrix with vector b and stores in dst
    fn matrix_vector_product(dst: &mut Vector2D, b: &Vector2D, a_plus_x: &Vector2D, a_plus_y: &Vector2D, a_diag: &Vector2D) {
        for y in 0..dst.rows {
            for x in 0..dst.columns {
                let mut t = a_diag.at(y, x) * b.at(y, x);

                if x > 0 {
                    t += a_plus_x.at(y, x - 1) * b.at(y, x - 1);
                }
                if y > 0 {
                    t += a_plus_y.at(y - 1, x) * b.at(y - 1, x);
                }
                if x < (dst.columns - 1) {
                    t += a_plus_x.at(y, x) * b.at(y, x + 1);
                }
                if y < (dst.rows - 1) {
                    t += a_plus_y.at(y, x) * b.at(y + 1, x);
                }

                *dst.at_mut(y, x) = t;
            }
        }
    }

    fn scaled_add1(dst: &mut Vector2D, b: &Vector2D, s: f64) {
        for y in 0..dst.rows {
            for x in 0..dst.columns {
                *dst.at_mut(y, x) += b.at(y, x) * s;
            }
        }
    }

    fn scaled_add2(dst: &mut Vector2D, a: &Vector2D,  s: f64) {
        for y in 0..dst.rows {
            for x in 0..dst.columns {
                *dst.at_mut(y, x) = a.at(y, x) + dst.at(y, x) * s;
            }
        }
    }

    fn infinity_norm(a: &Vector2D) -> f64 {
        let mut max_a = 0.0;
        for y in 0..a.rows {
            for x in 0..a.columns {
                max_a = max(max_a, (a.at(y, x)).abs());
            }
        }
        max_a
    }

    // Solves and mutates pressure array based on divergence, passed to linear solver function
//    fn solve_pressure(&mut self) {
//        (self.linear_solver)(&mut self.pressure, &self.r, self.fluid_density, self.dt, self.dx, 100);
//    }

    // Applies computed pressure field to the xy velocity vector field
    fn apply_pressure(&mut self) {
        let scale = self.dt / (self.fluid_density * self.dx);

        for row in 0..self.rows {
            for column in 0..self.columns {
                *self.x_velocity_src.at_mut(row, column) -= scale * self.p.at(row, column);
                *self.x_velocity_src.at_mut(row, column + 1) += scale * self.p.at(row, column);
                *self.y_velocity_src.at_mut(row, column) -= scale * self.p.at(row, column);
                *self.y_velocity_src.at_mut(row + 1, column) += scale * self.p.at(row, column);
            }
        }
    }

    // Sets boundaries of simulation by setting xy velocities at boundaries to 0
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


    // Projection method implements each step of the calculation
    fn project(&mut self) {
        self.set_boundaries();
        self.calculate_rhs();
        self.build_pressure_matrix();
        self.build_preconditioner();

        // Set initial guess of zeroes
        self.p.field = vec![0.0; self.rows * self.columns];


        // Apply preconditioner to z
        FluidSolver::apply_preconditioner(&mut self.z, &self.r, &self.a_plus_x, &self.a_plus_y, &self.pre_con);

        self.s.field = (self.z.field).clone();

        let mut max_error = FluidSolver::infinity_norm(&self.r);

        if max_error < 1e-5 {
            return;
        }

        let mut sigma = FluidSolver::dot_product(&self.z, &self.r);

        for iteration in 0..600 {
            FluidSolver::matrix_vector_product(&mut self.z, &self.s, &self.a_plus_x, &self.a_plus_y, &self.a_diag);

            let alpha = sigma / FluidSolver::dot_product(&self.z, &self.s);

            FluidSolver::scaled_add1(&mut self.p, &self.s, alpha);
            FluidSolver::scaled_add1(&mut self.r, &self.z, -alpha);

            max_error = FluidSolver::infinity_norm(&self.r);

            if max_error < 1e-5 {
                println!("Exiting solver after {} iterations, maximum error is {}", iteration, max_error);
                break;
            }

            FluidSolver::apply_preconditioner(&mut self.z, &self.r, &self.a_plus_x, &self.a_plus_y, &self.pre_con);

            let sigma_new = FluidSolver::dot_product(&self.z, &self.r);
            FluidSolver::scaled_add2(&mut self.s, &self.z, sigma_new / sigma);
            sigma = sigma_new;

            if iteration == 599 {
                println!("Exceeded budget of {} iterations, maximum change was {}", 600, max_error);
            }
        }

        self.apply_pressure();

        self.set_boundaries();
    }

    // Advection method moves density scalar field through velocity vector field to produce output
    fn advect(&mut self) {
        // TODO need to set up two different interpolation methods
        let interpolation = interpolation::bilinear_interpolate;
        (self.advection)(&mut self.x_velocity_dst, &self.x_velocity_src, &self.x_velocity_src, &self.y_velocity_src, self.dt, self.dx, &interpolation, &self.integration);
        (self.advection)(&mut self.y_velocity_dst, &self.y_velocity_src, &self.x_velocity_src, &self.y_velocity_src, self.dt, self.dx, &interpolation, &self.integration);
        (self.advection)(&mut self.density_dst, &self.density_src, &self.x_velocity_src, &self.y_velocity_src, self.dt, self.dx, &interpolation, &self.integration);
        self.swap();
    }

    // Produces the next frame of the simulation by projecting then advecting
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

    // Basic function to convert density_src array into an image buffer
    pub fn to_image(&self, buffer: &mut Vec<u8>) {
        for i in 0..(self.rows * self.columns) {
            let shade: u8 = (self.density_src.field[i] * 255.0 / 3.0) as u8;

            buffer[i * 4 + 0] = shade;
            buffer[i * 4 + 1] = shade;
            buffer[i * 4 + 2] = shade;
            buffer[i * 4 + 3] = 0xFF; // Alpha is 255
        }
    }
}