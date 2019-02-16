use crate::field::Field;
use crate::types::Interpolation;
use crate::types::Integration;

pub fn empty(_: &mut Field, _: &Field, _: &Field, _: f64, _: f64, _: &Interpolation, _: &Integration) {

}

pub fn semi_lagrangian(field: &mut Field, x_velocity: &Field, y_velocity: &Field, dt: f64, dx: f64, interpolation: &Interpolation, integration: &Integration) {
    for row in 0..field.rows {
        for column in 0..field.columns {
            let mut x = column as f64 + field.offset_x;
            let mut y = row as f64 + field.offset_y;

            let f1 = |_: f64, v: f64| (-interpolation(v, y, &x_velocity) / dx);
            let f2 = |_: f64, v: f64| (-interpolation(x, v, &y_velocity) / dx);

            let old_x = integration(0.0, x, &f1, dt);
            let old_y = integration(0.0, y, &f2, dt);

            *field.at_mut(row, column) = interpolation(old_x, old_y, field);
        }
    }
}