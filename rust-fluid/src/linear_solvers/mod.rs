use crate::field::Field;

pub fn empty(_: &mut Field, _: &Field, _: f64, _: f64, _: f64, _: usize) {

}

pub fn gauss_siedel(pressure: &mut Field, divergence: &Field, fluid_density: f64, dt: f64, dx: f64, limit: usize) {
    let rows = pressure.rows;
    let columns = pressure.columns;
    let scale = dt / (fluid_density * dx * dx);

    let mut max_delta = 0.0;

    for iteration in 0..limit {
        max_delta = 0.0;

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

//                println!("scale: {}, old: {}, new {}", scale, pressure.at(row, column), new_pressure);

                let delta = (pressure.at(row, column) - new_pressure).abs();



                if max_delta < delta {
                    max_delta = delta;
                }

                *pressure.at_mut(row, column) = new_pressure;
            }
        }
//        println!("{}", max_delta);
        if max_delta < 1e-5 {
            println!("Exiting solver after {} iterations, maximum change is {}", iteration, max_delta);
            return;
        }
    }
    println!("Exceeded budget of {} iterations, maximum change was {}", limit, max_delta);
}

pub fn relaxation(x: &mut Field, b: &Field, density: f64, dt: f64, dx: f64, limit: usize) {
    let w = x.columns;
    let h = x.rows;

    for _ in 0..limit {
        for r in 0..h {
            for c in 0..w {

                //     |p3|
                //  ---|--|---
                //  p1 |  | p2
                //  ---|--|---
                //     |p4|


                let mut alpha = 4.0;


                let p1 = if c as i32 - 1 >= 0 { x.at(r, c-1) } else { alpha -= 1.0; 0.0 } * (dt / ( density * dx * dx ));
                let p2 = if c as i32 + 1 < w as i32 { x.at(r, c+1) } else { alpha-=1.0; 0.0 } * (dt / ( density * dx * dx ));
                let p3 = if r as i32 + 1 < h as i32 { x.at(r+1, c) } else { alpha-=1.0; 0.0 } * (dt / ( density * dx * dx ));
                let p4 = if r as i32 - 1 >= 0 { x.at(r-1, c) } else { alpha-=1.0; 0.0 } * (dt / ( density * dx * dx ));

                *x.at_mut(r, c) = ( b.at(r, c) + p1 + p2 + p3 + p4 ) / (alpha * (dt / ( density * dx * dx )));
            }
        }
    }
}