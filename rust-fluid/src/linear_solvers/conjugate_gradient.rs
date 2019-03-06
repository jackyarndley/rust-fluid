use crate::util::{Sparse, infinity_norm, dot_product, matrix_vector_product, scaled_add1, scaled_add2};

fn build_pressure_matrix(a: &mut Sparse, fluid_density: f64, dt: f64, dx: f64, rows: usize, columns: usize) {
    let scale = dt / (fluid_density * dx * dx);
    a.diagonals = vec![0.0; rows * columns];

    for row in 0..rows {
        for column in 0..columns {
            let element = row * columns + column;

            if column < columns - 1 {
                a.diagonals[element] += scale;
                a.diagonals[element + 1] += scale;
                a.plus_x[element] = -scale;
            } else {
                a.plus_x[element] = 0.0;
            }

            if row < rows - 1 {
                a.diagonals[element] += scale;
                a.diagonals[element + columns] += scale;
                a.plus_y[element] = -scale;
            } else {
                a.plus_y[element] = 0.0;
            }
        }
    }
}

fn build_preconditioner(preconditioner: &mut Vec<f64>, a: &Sparse, rows: usize, columns: usize) {
    let tau = 0.97;

    for row in 0..rows {
        for column in 0..columns {
            let element = row * columns + column;

            let mut e = a.diagonals[element];

            if column > 0 {
                let px = a.plus_x[element - 1] * preconditioner[element - 1];
                let py = a.plus_y[element - 1] * preconditioner[element - 1];
                e -= px * px + tau * px * py;
            }

            if row > 0 {
                let px = a.plus_x[element - columns] * preconditioner[element - columns];
                let py = a.plus_y[element - columns] * preconditioner[element - columns];
                e -= py * py + tau * px * py;
            }

            preconditioner[element] = 1.0 / (e + 1e-30).sqrt();
        }
    }
}

fn apply_preconditioner(auxiliary: &mut Vec<f64>, residual: &Vec<f64>, a: &Sparse, preconditioner: &Vec<f64>, rows: usize, columns: usize) {
    for row in 0..rows {
        for column in 0..columns {
            let element = row * columns + column;
            let mut t = residual[element];

            if column > 0 {
                t -= a.plus_x[element - 1] * preconditioner[element - 1] * auxiliary[element - 1];
            }

            if row > 0 {
                t -= a.plus_y[element - columns] * preconditioner[element - columns] * auxiliary[element - columns];
            }

            auxiliary[element] = t * preconditioner[element];
        }
    }

    for row in (0..rows).rev() {
        for column in (0..columns).rev() {
            let element = row * columns + column;
            let mut t = auxiliary[element];

            if column < columns - 1 {
                t -= a.plus_x[element] * preconditioner[element] * auxiliary[element + 1];
            }

            if row < rows - 1 {
                t -= a.plus_y[element] * preconditioner[element] * auxiliary[element + columns];
            }

            auxiliary[element] = t * preconditioner[element];
        }
    }
}

pub fn conjugate_gradient(pressure: &mut Vec<f64>,
                          residual: &mut Vec<f64>,
                          auxiliary: &mut Vec<f64>,
                          search: &mut Vec<f64>,
                          preconditioner: &mut Vec<f64>,
                          a: &mut Sparse,
                          fluid_density: f64,
                          dt: f64,
                          dx: f64,
                          rows: usize,
                          columns: usize,
                          limit: usize) {

    build_pressure_matrix(a, fluid_density, dt, dx, rows, columns);
    build_preconditioner(preconditioner, a, rows, columns);

    *pressure = vec![0.0; rows * columns];

    apply_preconditioner(auxiliary, residual, a, preconditioner, rows, columns);
    *search = auxiliary.clone();

    let mut max_error = infinity_norm(residual);

    if max_error < 1e-5 {
        return;
    }

    let mut sigma = dot_product(auxiliary, residual);

    for iteration in 0..limit {
        matrix_vector_product(auxiliary, search, a, rows, columns);

        let alpha = sigma / dot_product(auxiliary, search);

        scaled_add1(pressure, search, alpha);
        scaled_add1(residual, auxiliary, -alpha);

        max_error = infinity_norm(residual);

        if max_error < 1e-5 {
            println!("Exiting solver after {} iterations, maximum change is {}", iteration, max_error);
            return;
        }

        apply_preconditioner(auxiliary, residual, a, preconditioner, rows, columns);

        let sigma_new = dot_product(auxiliary, residual);
        scaled_add2(search, auxiliary, sigma_new / sigma);
        sigma = sigma_new;
    }
    println!("Exceeded budget of {} iterations, maximum change was {}", limit, max_error)
}