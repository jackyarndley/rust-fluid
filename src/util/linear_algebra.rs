use crate::util::sparse::Sparse;
use crate::util::helper::max;

pub fn dot_product(a: &Vec<f64>, b: &Vec<f64>) -> f64 {
    let mut result = 0.0;
    for element in 0..a.len() {
        result += a[element] * b[element];
    }
    result
}

// Multiplies pressure matrix with vector b and stores in dst
pub fn matrix_vector_product(dst: &mut Vec<f64>, b: &Vec<f64>, a: &Sparse, rows: usize, columns: usize) {
    for row in 0..rows {
        for column in 0..columns {
            let element = row * columns + column;
            let mut t = a.diagonals[element] * b[element];

            if column > 0 {
                t += a.plus_x[element - 1] * b[element - 1];
            }

            if row > 0 {
                t += a.plus_y[element - columns] * b[element - columns];
            }

            if column < columns - 1 {
                t += a.plus_x[element] * b[element + 1];
            }

            if row < rows - 1 {
                t += a.plus_y[element] * b[element + columns];
            }

            dst[element] = t;
        }
    }
}

pub fn scaled_add1(dst: &mut Vec<f64>, b: &Vec<f64>, s: f64) {
    for element in 0..dst.len() {
        dst[element] += b[element] * s;
    }
}

pub fn scaled_add2(dst: &mut Vec<f64>, a: &Vec<f64>,  s: f64) {
    for element in 0..dst.len() {
        dst[element] = a[element] + dst[element] * s;
    }
}

pub fn infinity_norm(a: &Vec<f64>) -> f64 {
    let mut max_a = 0.0;
    for element in a {
        max_a = max(max_a, element.abs());
    }
    max_a
}