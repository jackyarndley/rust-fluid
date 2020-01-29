use crate::util::sparse::Sparse;
use crate::util::linear_algebra::{infinity_norm, dot_product, matrix_vector_product, scaled_add1, scaled_add2};
use crate::util::fluid_quantity::FluidQuantity;

fn build_pressure_matrix(a: &mut Sparse, cell: &Vec<u8>, fluid_density: f64, dt: f64, dx: f64, rows: usize, columns: usize, u_velocity: &FluidQuantity, v_velocity: &FluidQuantity) {
    let scale = dt / (fluid_density * dx * dx);
    a.diagonals = vec![0.0; rows * columns];
    a.plus_x = vec![0.0; rows * columns];
    a.plus_y = vec![0.0; rows * columns];

    for row in 0..rows {
        for column in 0..columns {
            let element = row * columns + column;
            if cell[element] == 0 {
                if column < columns - 1 && cell[element + 1] == 0 {
                    let factor = scale * u_velocity.volume_at(row, column + 1);
                    a.diagonals[element] += factor;
                    a.diagonals[element + 1] += factor;
                    a.plus_x[element] = -factor;
                }

                if row < rows - 1 && cell[element + columns] == 0 {
                    let factor = scale * v_velocity.volume_at(row + 1, column);
                    a.diagonals[element] += factor;
                    a.diagonals[element + columns] += factor;
                    a.plus_y[element] = -factor;
                }
            }
        }
    }
}

fn build_preconditioner(preconditioner: &mut Vec<f64>, a: &Sparse, cell: &Vec<u8>, rows: usize, columns: usize) {
    let tau = 0.97;

    for row in 0..rows {
        for column in 0..columns {
            let element = row * columns + column;

            if cell[element] == 0 {
                let mut e = a.diagonals[element];

                if column > 0 && cell[element - 1] == 0 {
                    let px = a.plus_x[element - 1] * preconditioner[element - 1];
                    let py = a.plus_y[element - 1] * preconditioner[element - 1];
                    e -= px * px + tau * px * py;
                }

                if row > 0 && cell[element - columns] == 0 {
                    let px = a.plus_x[element - columns] * preconditioner[element - columns];
                    let py = a.plus_y[element - columns] * preconditioner[element - columns];
                    e -= py * py + tau * px * py;
                }

                preconditioner[element] = 1.0 / (e + 1e-30).sqrt();
            }
        }
    }
}

fn apply_preconditioner(auxiliary: &mut Vec<f64>, residual: &Vec<f64>, a: &Sparse, preconditioner: &Vec<f64>, cell: &Vec<u8>, rows: usize, columns: usize) {
    for row in 0..rows {
        for column in 0..columns {
            let element = row * columns + column;

            if cell[element] == 0 {
                let mut t = residual[element];

                if column > 0 && cell[element - 1] == 0 {
                    t -= a.plus_x[element - 1] * preconditioner[element - 1] * auxiliary[element - 1];
                }

                if row > 0 && cell[element - columns] == 0 {
                    t -= a.plus_y[element - columns] * preconditioner[element - columns] * auxiliary[element - columns];
                }

                auxiliary[element] = t * preconditioner[element];
            }
        }
    }

    for row in (0..rows).rev() {
        for column in (0..columns).rev() {
            let element = row * columns + column;

            if cell[element] == 0 {
                let mut t = auxiliary[element];

                if column < columns - 1 && cell[element + 1] == 0 {
                    t -= a.plus_x[element] * preconditioner[element] * auxiliary[element + 1];
                }

                if row < rows - 1 && cell[element + columns] == 0 {
                    t -= a.plus_y[element] * preconditioner[element] * auxiliary[element + columns];
                }

                auxiliary[element] = t * preconditioner[element];
            }
        }
    }
}

pub fn conjugate_gradient(pressure: &mut Vec<f64>,
                          residual: &mut Vec<f64>,
                          auxiliary: &mut Vec<f64>,
                          search: &mut Vec<f64>,
                          preconditioner: &mut Vec<f64>,
                          a: &mut Sparse,
                          cell: &Vec<u8>,
                          fluid_density: f64,
                          dt: f64,
                          dx: f64,
                          rows: usize,
                          columns: usize,
                          limit: usize,
                          u_velocity: &FluidQuantity,
                          v_velocity: &FluidQuantity) {

    build_pressure_matrix(a, cell, fluid_density, dt, dx, rows, columns, u_velocity, v_velocity);
    build_preconditioner(preconditioner, a, cell, rows, columns);

    *pressure = vec![0.0; rows * columns];

    apply_preconditioner(auxiliary, residual, a, preconditioner, cell, rows, columns);
    *search = auxiliary.clone();

    let mut max_error = infinity_norm(residual);

    if max_error < 1e-4 {
        return;
    }

    let mut sigma = dot_product(auxiliary, residual);

    for _iteration in 0..limit {
        matrix_vector_product(auxiliary, search, a, rows, columns);

        let alpha = sigma / dot_product(auxiliary, search);

        scaled_add1(pressure, search, alpha);
        scaled_add1(residual, auxiliary, -alpha);

        max_error = infinity_norm(residual);

        if max_error < 1e-4 {
            return;
        }

        apply_preconditioner(auxiliary, residual, a, preconditioner, cell, rows, columns);

        let sigma_new = dot_product(auxiliary, residual);
        scaled_add2(search, auxiliary, sigma_new / sigma);
        sigma = sigma_new;
    }
    println!("Solver exceeded budget of {} iterations, maximum change was {}", limit, max_error)
}