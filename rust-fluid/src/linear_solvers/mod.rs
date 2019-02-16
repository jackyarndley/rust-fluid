use crate::field::Field;

pub fn empty(_: &mut Field, _: &Field, _: f64, _: f64, _: f64, _: usize) {

}

pub fn gauss_siedel(pressure: &mut Field, divergence: &Field, fluid_density: f64, dt: f64, dx: f64, limit: usize) {
    let rows = pressure.rows;
    let columns = pressure.columns;
    let scale = dt / (fluid_density * dx * dx);

    for _iteration in 0..limit {
        let mut max_delta = 0.0;

        for row in 0..rows {
            for column in 0..columns {
                let mut diag = 0.0;
                let mut off_diag = 0.0;

                if column > 0 {
                    diag += scale;
                    off_diag -= scale * pressure.at(row, column - 1);
                }

                if row > 0 {
                    diag += scale;
                    off_diag -= scale * pressure.at(row - 1, column);
                }

                if column < (columns - 1) {
                    diag += scale;
                    off_diag -= scale * pressure.at(row, column + 1);
                }

                if row < (rows - 1) {
                    diag += scale;
                    off_diag -= scale * pressure.at(row + 1, column);
                }

                let new_pressure = (divergence.at(row, column) - off_diag) / diag;

                let delta = (pressure.at(row, column) - new_pressure).abs();

                if max_delta < delta {
                    max_delta = delta;
                }

                *pressure.at_mut(row, column) = new_pressure;
            }
        }
    }
}