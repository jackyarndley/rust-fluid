pub fn gauss_siedel(pressure: &mut Vec<f64>, residual: &mut Vec<f64>, fluid_density: f64, dt: f64, dx: f64, rows: usize, columns: usize, limit: usize) {
    let scale = dt / (fluid_density * dx * dx);
    let mut max_delta = 0.0;

    for iteration in 0..limit {
        max_delta = 0.0;

        for row in 0..rows {
            for column in 0..columns {
                let mut diagonal = 0.0;
                let mut off_diagonal = 0.0;
                let element = row * columns + column;

                if column > 0 {
                    diagonal += scale;
                    off_diagonal -= scale * pressure[element - 1];
                }

                if row > 0 {
                    diagonal += scale;
                    off_diagonal -= scale * pressure[element - columns];
                }

                if column < (columns - 1) {
                    diagonal += scale;
                    off_diagonal -= scale * pressure[element + 1];
                }

                if row < (rows - 1) {
                    diagonal += scale;
                    off_diagonal -= scale * pressure[element + columns];
                }

                let new_pressure = (residual[element] - off_diagonal) / diagonal;
                let delta = (pressure[element] - new_pressure).abs();
                if max_delta < delta {
                    max_delta = delta;
                };

                pressure[element] = new_pressure;
            }
        }

        if max_delta < 1e-5 {
            println!("Exiting solver after {} iterations, maximum change is {}", iteration, max_delta);
            return;
        }
    }
    println!("Exceeded budget of {} iterations, maximum change was {}", limit, max_delta);
}