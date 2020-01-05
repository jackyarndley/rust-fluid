use crate::util::Field;
use crate::interpolation::bilinear_interpolate;
use crate::integration::Integration;

//pub fn semi_lagrangian(field_dst: &mut Field, field_src: &Field, x_velocity: &Field, y_velocity: &Field, dt: f64, dx: f64, ) {
////    for row in 0..field_dst.rows {
////        for column in 0..field_dst.columns {
////            let x = column as f64 + field_dst.x_offset;
////            let y = row as f64 + field_dst.y_offset;
////
////            let f1 = |_t: f64, y0: f64| -1.0 * bilinear_interpolate(y0, y, &x_velocity) / dx;
////            let f2 = |_t: f64, y0: f64| -1.0 * bilinear_interpolate(x, y0, &y_velocity) / dx;
////
////            let old_x: f64 = integration(0.0, x, &f1, dt);
////            let old_y: f64 = integration(0.0, y, &f2, dt);
////
////            *field_dst.at_mut(row, column) = bilinear_interpolate(old_x, old_y, field_src);
////        }
////    }
//}

pub enum Advection {
    Empty,
    SemiLagrangian
}

impl Advection {
    pub fn advect(&self, field_dst: &mut Field, field_src: &Field, x_velocity: &Field, y_velocity: &Field, dt: f64, dx: f64, integration: &Integration) {
        match self {
            Advection::Empty => {

            }
            Advection::SemiLagrangian => {
                for row in 0..field_dst.rows {
                    for column in 0..field_dst.columns {
                        let x = column as f64 + field_dst.x_offset;
                        let y = row as f64 + field_dst.y_offset;

                        let f1 = |_t: f64, y0: f64| -1.0 * bilinear_interpolate(y0, y, &x_velocity) / dx;
                        let f2 = |_t: f64, y0: f64| -1.0 * bilinear_interpolate(x, y0, &y_velocity) / dx;

                        let old_x: f64 = integration.run(0.0, x, &f1, dt);
                        let old_y: f64 = integration.run(0.0, y, &f2, dt);

                        *field_dst.at_mut(row, column) = bilinear_interpolate(old_x, old_y, field_src);
                    }
                }
            }
        }
    }
}