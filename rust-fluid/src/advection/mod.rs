use crate::field::Field;
use crate::types::Interpolation;
use crate::types::Integration;

pub fn empty(_: &mut Field, _: &Field, _: &Field, _: f64, _: f64, _: &Fn(f64, f64, &Field) -> f64, _: &Fn(f64, f64, &Fn(f64, f64) -> f64, f64) -> f64) {

}

pub fn semi_lagrangian(field: &mut Field, x_velocity: &Field, y_velocity: &Field, dt: f64, dx: f64, interpolation: &Fn(f64, f64, &Field) -> f64, integration: &Fn(f64, f64, &Fn(f64, f64) -> f64, f64) -> f64) {
    for row in 0..field.rows {
        for column in 0..field.columns {
            let x = column as f64 + field.offset_x;
            let y = row as f64 + field.offset_y;

            let f1 = |_: f64, o: f64| -interpolation(o, y, &x_velocity) / dx;
            let f2 = |_: f64, o: f64| -interpolation(x, o, &y_velocity) / dx;

            let old_x: f64 = integration(0.0, x, &f1, dt);
            let old_y: f64 = integration(0.0, y, &f2, dt);

            *field.at_mut(row, column) = interpolation(old_x, old_y, field);
        }
    }
}