use crate::grid::Grid;
use crate::functions::*;
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

    pub fn add_source2(&mut self, row1: usize, column1: usize, row2: usize, column2: usize, velocity_x: f64, velocity_y: f64, density: f64) {
        let f_row1 = (row1 as f64) / (self.rows as f64);
        let f_row2 = (row2 as f64) / (self.rows as f64);
        let f_column1 = (column1 as f64) / (self.columns as f64);
        let f_column2 = (column2 as f64) / (self.columns as f64);

        // println!("{} {} {} {}", f_row1, f_row2, f_column1, f_column2);

        for row in row1..row2 {
            for column in column1..column2 {

                let length = length(
                    (2.0 * (column as f64 + 0.5) * self.dx - (f_column1 + f_column2)) / (f_column2 - f_column1),
                    (2.0 * (row as f64 + 0.5) * self.dx - (f_row1 + f_row2)) / (f_row2 - f_row1)
                );
                // println!("length {}", length);

                let vi_velocity_x = cubic_pulse(length) * velocity_x;
                let vi_velocity_y = cubic_pulse(length) * velocity_y;
                let vi_density = cubic_pulse(length) * density;

                if self.x_velocity_src.at(row, column).abs() < vi_velocity_x.abs() {
                    *self.x_velocity_src.at_mut(row, column) = vi_velocity_x;
                }

                if self.y_velocity_src.at(row, column).abs() < vi_velocity_y.abs() {
                    *self.y_velocity_src.at_mut(row, column) = vi_velocity_y;
                }
                // println!("density {}", vi_density);
                if self.density_src.at(row, column).abs() < vi_density.abs() {
                    // println!("density {}", vi_density);
                    *self.density_src.at_mut(row, column) = vi_density;
                }
            }
        }
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

                self.runge_kutta3(&mut density_x, &mut density_y, time_step);

                *self.density_dst.at_mut(row, column) = cerp_grid(&self.density_src, density_x, density_y);
            }
        }

        // Advect x_velocity Grid
        for row in 0..self.x_velocity_src.rows {
            for column in 0..self.x_velocity_src.columns {
                let mut x_velocity_x = column as f64 + self.x_velocity_src.offset_x;
                let mut x_velocity_y = row as f64 + self.x_velocity_src.offset_y;

                self.runge_kutta3(&mut x_velocity_x, &mut x_velocity_y, time_step);

                *self.x_velocity_dst.at_mut(row, column) = cerp_grid(&self.x_velocity_src, x_velocity_x, x_velocity_y);
            }
        }

        // Advect y_velocity Grid
        for row in 0..self.y_velocity_src.rows {
            for column in 0..self.y_velocity_src.columns {
                let mut y_velocity_x = column as f64 + self.y_velocity_src.offset_x;
                let mut y_velocity_y = row as f64 + self.y_velocity_src.offset_y;;

                self.runge_kutta3(&mut y_velocity_x, &mut y_velocity_y, time_step);

                *self.y_velocity_dst.at_mut(row, column) = cerp_grid(&self.y_velocity_src, y_velocity_x, y_velocity_y);
            }
        }
    }

//    fn euler(&mut self, x: &mut f64, y: &mut f64, time_step: f64) {
//        let x_velocity = lerp_grid(&self.x_velocity_src, *x, *y) / self.dx;
//        let y_velocity = lerp_grid(&self.y_velocity_src, *x, *y) / self.dx;
//
//        *x -= x_velocity * time_step;
//        *y -= y_velocity * time_step;
//    }

    fn runge_kutta3(&mut self, x: &mut f64, y: &mut f64, time_step: f64) {
        let x_velocity1 = lerp_grid(&self.x_velocity_src, *x, *y) / self.dx;
        let y_velocity1 = lerp_grid(&self.y_velocity_src, *x, *y) / self.dx;

        let mid_x = *x - 0.5 * time_step * x_velocity1;
        let mid_y = *y - 0.5 * time_step * y_velocity1;

        let x_velocity2 = lerp_grid(&self.x_velocity_src, mid_x, mid_y) / self.dx;
        let y_velocity2 = lerp_grid(&self.y_velocity_src, mid_x, mid_y) / self.dx;

        let last_x = *x - 0.75 * time_step * x_velocity2;
        let last_y = *y - 0.75 * time_step * y_velocity2;

        let x_velocity3 = lerp_grid(&self.x_velocity_src, last_x, last_y) / self.dx;
        let y_velocity3 = lerp_grid(&self.y_velocity_src, last_x, last_y) / self.dx;

        *x -= time_step * ((2.0/9.0) * x_velocity1 + (3.0/9.0) * x_velocity2 + (4.0/9.0) * x_velocity3);
        *y -= time_step * ((2.0/9.0) * y_velocity1 + (3.0/9.0) * y_velocity2 + (4.0/9.0) * y_velocity3);
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