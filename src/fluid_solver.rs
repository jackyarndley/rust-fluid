use crate::grid::Grid;
use std::mem::swap;

pub struct FluidSolver {
    pub x_velocity_src: Grid,
    pub x_velocity_dst: Grid,
    pub y_velocity_src: Grid,
    pub y_velocity_dst: Grid,
    pub density_src: Grid,
    pub density_dst: Grid,
    pub divergence: Grid,
    pub pressure: Grid,
    pub rows: usize,
    pub columns: usize,
    pub fluid_density: f64,
    pub dx: f64
}

impl FluidSolver {
    // Creates a new FluidSolver struct with parameters supplied
    pub fn new(rows: usize, columns: usize, fluid_density: f64) -> FluidSolver {
        FluidSolver {
            x_velocity_src: Grid::new(rows, columns + 1, 0.0, 0.5),
            x_velocity_dst: Grid::new(rows, columns + 1, 0.0, 0.5),
            y_velocity_src: Grid::new(rows + 1, columns, 0.5, 0.0),
            y_velocity_dst: Grid::new(rows + 1, columns, 0.5, 0.0),
            density_src: Grid::new(rows, columns, 0.5, 0.5),
            density_dst: Grid::new(rows, columns, 0.5, 0.5),
            divergence: Grid::new(rows, columns, 0.5, 0.5),
            pressure: Grid::new(rows, columns, 0.5, 0.5),
            rows,
            columns,
            fluid_density,
            dx: 1.0 / rows as f64
        }
    }

    // Sets an individual cell to have specified velocity and density
    pub fn add_source(&mut self, row: usize, column: usize, velocity_x: f64, velocity_y: f64, density: f64) {
        *self.x_velocity_src.at_mut(row, column) = velocity_x;
        *self.y_velocity_src.at_mut(row, column) = velocity_y;
        *self.density_src.at_mut(row, column) = density;
    }

    pub fn update(&mut self, time_step: f64) {
        // Calculate divergence of new velocity field
        self.calculate_divergence();

        // Solve for new pressure field
        self.solve_pressure(100, time_step);

        // Apply pressure field to velocity field
        self.apply_pressure(time_step);

        // Advect each Grid
        self.advect(time_step);

        // Swap buffers
        self.swap_grid();
    }

    // Calculates divergence of velocity field
    fn calculate_divergence(&mut self) {
        let scale: f64 = 1.0 / self.dx;

        for row in 0..self.rows {
            for column in 0..self.columns {
                let x_velocity1 = self.x_velocity_src.at(row, column);
                let x_velocity2 = self.x_velocity_src.at(row, column + 1);
                let y_velocity1 = self.y_velocity_src.at(row, column);
                let y_velocity2 = self.y_velocity_src.at(row + 1, column);

                *self.divergence.at_mut(row, column) = -1.0 * scale * (x_velocity2 - x_velocity1 + y_velocity2 - y_velocity1);
            }
        }
    }

    // Pressure solve using 'Gauss-Siedel'
    fn solve_pressure(&mut self, limit: usize, time_step: f64) {
        let scale: f64 = time_step / (self.fluid_density * self.dx * self.dx);
        let mut max_delta = 0.0;

        for iteration in 0..limit {
            max_delta = 0.0;

            for row in 0..self.rows {
                for column in 0..self.columns {
                    let mut diag = 0.0;
                    let mut off_diag = 0.0;

                    if column > 0 {
                        diag += scale;
                        off_diag -= scale * self.pressure.at(row, column - 1);
                    }

                    if row > 0 {
                        diag += scale;
                        off_diag -= scale * self.pressure.at(row - 1, column);
                    }

                    if column < (self.columns - 1) {
                        diag += scale;
                        off_diag -= scale * self.pressure.at(row, column + 1);
                    }

                    if row < (self.rows - 1) {
                        diag += scale;
                        off_diag -= scale * self.pressure.at(row + 1, column);
                    }

                    let new_pressure = (self.divergence.at(row, column) - off_diag) / diag;

                    let delta = (self.pressure.at(row, column) - new_pressure).abs();

                    if max_delta < delta {
                        max_delta = delta;
                    }

                    *self.pressure.at_mut(row, column) = new_pressure;
                }
            }

            if max_delta < 1e-5 {

                println!("Exiting solver after {} iterations, maximum change is {}", iteration, max_delta);
                return;
            }
        }

        println!("Exceeded budget of {} iterations, maximum change was {}", limit, max_delta);
    }

    // Applies computed pressure field to velocity field
    fn apply_pressure(&mut self, time_step: f64) {
        let scale: f64 = time_step / (self.fluid_density * self.dx);

        for row in 0..self.rows {
            for column in 0..self.columns {
                *self.x_velocity_src.at_mut(row, column) -= scale * self.pressure.at(row, column);
                *self.x_velocity_src.at_mut(row, column + 1) += scale * self.pressure.at(row, column);
                *self.y_velocity_src.at_mut(row, column) -= scale * self.pressure.at(row, column);
                *self.y_velocity_src.at_mut(row + 1, column) += scale * self.pressure.at(row, column);
            }
        }

        for row in 0..self.rows {
            *self.x_velocity_src.at_mut(row, 0) = 0.0;
            *self.x_velocity_src.at_mut(row, self.columns) = 0.0;
        }

        for column in 0..self.columns {
            *self.y_velocity_src.at_mut(0, column) = 0.0;
            *self.y_velocity_src.at_mut(self.rows, column) = 0.0;
        }
    }

    fn advect(&mut self, time_step: f64) {
        // Advect density Grid
        for row in 0..self.density_src.rows {
            for column in 0..self.density_src.columns {
                let mut density_x = column as f64 + self.density_src.offset_x;
                let mut density_y = row as f64 + self.density_src.offset_y;

                self.euler(&mut density_x, &mut density_y, time_step);

                *self.density_dst.at_mut(row, column) = FluidSolver::lerp_grid(&self.density_src, density_x, density_y);
            }
        }

        // Advect x_velocity Grid
        for row in 0..self.x_velocity_src.rows {
            for column in 0..self.x_velocity_src.columns {
                let mut x_velocity_x = column as f64 + self.x_velocity_src.offset_x;
                let mut x_velocity_y = row as f64 + self.x_velocity_src.offset_y;

                self.euler(&mut x_velocity_x, &mut x_velocity_y, time_step);

                *self.x_velocity_dst.at_mut(row, column) = FluidSolver::lerp_grid(&self.x_velocity_src, x_velocity_x, x_velocity_y);
            }
        }

        // Advect y_velocity Grid
        for row in 0..self.y_velocity_src.rows {
            for column in 0..self.y_velocity_src.columns {
                let mut y_velocity_x = column as f64 + self.y_velocity_src.offset_x;
                let mut y_velocity_y = row as f64 + self.y_velocity_src.offset_y;;

                self.euler(&mut y_velocity_x, &mut y_velocity_y, time_step);

                *self.y_velocity_dst.at_mut(row, column) = FluidSolver::lerp_grid(&self.y_velocity_src, y_velocity_x, y_velocity_y);
            }
        }
    }

    fn euler(&mut self, x: &mut f64, y: &mut f64, time_step: f64) {
        let x_velocity = FluidSolver::lerp_grid(&self.x_velocity_src, *x, *y) / self.dx;
        let y_velocity = FluidSolver::lerp_grid(&self.y_velocity_src, *x, *y) / self.dx;

        *x -= x_velocity * time_step;
        *y -= y_velocity * time_step;
    }

    fn lerp_grid(grid: &Grid, x: f64, y: f64) -> f64 {
        let mut temp_x = x - grid.offset_x;
        let mut temp_y = y - grid.offset_y;

        if temp_x < 0.0 {
            temp_x = 0.0;
        } else if temp_x > (grid.columns as f64 - 1.001) {
            temp_x = grid.columns as f64 - 1.001;
        }

        if temp_y < 0.0 {
            temp_y = 0.0;
        } else if temp_y > (grid.rows as f64 - 1.001) {
            temp_y = grid.rows as f64 - 1.001;
        }

        let ix = temp_x as usize;
        let iy = temp_y as usize;

        temp_x -= ix as f64;
        temp_y -= iy as f64;

        let x00: f64 = grid.at(iy + 0, ix + 0);
        let x10: f64 = grid.at(iy + 0, ix + 1);
        let x01: f64 = grid.at(iy + 1, ix + 0);
        let x11: f64 = grid.at(iy + 1, ix + 1);

        FluidSolver::lerp(FluidSolver::lerp(x00, x10, temp_x), FluidSolver::lerp(x01, x11, temp_x), temp_y)
    }

    fn lerp(a: f64, b: f64, c: f64) -> f64 {
        a * (1.0 - c) + b * c
    }

    // Swaps src and dst grids
    fn swap_grid(&mut self) {
        swap(&mut self.density_src, &mut self.density_dst);
        swap(&mut self.x_velocity_src, &mut self.x_velocity_dst);
        swap(&mut self.y_velocity_src, &mut self.y_velocity_dst);
    }

    // Export Grid to image buffer
    pub fn to_image(&self, buffer: &mut Vec<u8>) {
        for i in 0..(self.rows * self.columns) {
            let shade: u8 = (255.0 * self.density_src.grid[i]) as u8;

            buffer[i * 4 + 0] = shade;
            buffer[i * 4 + 1] = shade;
            buffer[i * 4 + 2] = shade;
            buffer[i * 4 + 3] = 0xFF; // Alpha is 255
        }
    }
}