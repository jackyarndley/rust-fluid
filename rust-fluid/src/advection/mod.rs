use crate::field::Field;

pub fn empty(_: &mut Field, _: &Field, _: &Field, _: &Field, _: f64, _: f64, _: &Fn(f64, f64, &Field) -> f64, _: &Fn(f64, f64, &Fn(f64, f64) -> f64, f64) -> f64) {

}

pub fn semi_lagrangian(field_dst: &mut Field, field_src: &Field, x_velocity: &Field, y_velocity: &Field, dt: f64, dx: f64, interpolation: &Fn(f64, f64, &Field) -> f64, integration: &Fn(f64, f64, &Fn(f64, f64) -> f64, f64) -> f64) {
    for row in 0..field_dst.rows {
        for column in 0..field_dst.columns {
            let x = column as f64 + field_dst.offset_x;
            let y = row as f64 + field_dst.offset_y;

            let f1 = |_t: f64, y0: f64| -1.0 * interpolation(y0, y, &x_velocity) / dx;
            let f2 = |_t: f64, y0: f64| -1.0 * interpolation(x, y0, &y_velocity) / dx;

            let old_x: f64 = integration(0.0, x, &f1, dt);
            let old_y: f64 = integration(0.0, y, &f2, dt);

            *field_dst.at_mut(row, column) = interpolation(old_x, old_y, field_src);
        }
    }
}