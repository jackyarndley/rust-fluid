use crate::util::{Field, clamp, min, max, linear_interpolate, cubic_interpolate};

pub fn empty(_: f64, _: f64, _: &Field) -> f64 {
    0.0
}

// Bilinear interpolation takes 4 points surrounding a value and linearly interpolates their values
pub fn bilinear_interpolate(mut x: f64, mut y: f64, field: &Field) -> f64 {
    x = clamp(x - field.x_offset, 0.0, field.columns as f64 - 1.001);
    y = clamp(y - field.y_offset, 0.0, field.rows as f64 - 1.001);

    let p1_x = x.trunc() as usize;
    let p1_y = y.trunc() as usize;

    let p1 = field.at(p1_y, p1_x);
    let p2 = field.at(p1_y, p1_x + 1);
    let p3 = field.at(p1_y + 1, p1_x);
    let p4 = field.at(p1_y + 1, p1_x + 1);

    let l1 = linear_interpolate(p1, p2, x.fract());
    let l2 = linear_interpolate(p3, p4, x.fract());

    linear_interpolate(l1, l2, y.fract())
}

// Bicubic interpolation takes 16 points surrounding a value and cubic interpolates their values
pub fn bicubic_interpolate(mut x: f64, mut y: f64, field: &Field) -> f64 {
    x = clamp(x - field.x_offset, 0.0, field.columns as f64 - 1.001);
    y = clamp(y - field.y_offset, 0.0, field.rows as f64 - 1.001);

    let p1_x = x.trunc() as usize;
    let p1_y = y.trunc() as usize;

    let x0 = max(p1_x as isize - 1, 0) as usize;
    let x1 = p1_x;
    let x2 = p1_x + 1;
    let x3 = min(p1_x + 2, field.columns - 1);

    let y0 = max(p1_y as isize - 1, 0) as usize;
    let y1 = p1_y;
    let y2 = p1_y + 1;
    let y3 = min(p1_y + 2, field.rows - 1);

    let q0 = cubic_interpolate(field.at(y0, x0), field.at(y0, x1), field.at(y0, x2), field.at(y0, x3), x);
    let q1 = cubic_interpolate(field.at(y1, x0), field.at(y1, x1), field.at(y1, x2), field.at(y1, x3), x);
    let q2 = cubic_interpolate(field.at(y2, x0), field.at(y2, x1), field.at(y2, x2), field.at(y2, x3), x);
    let q3 = cubic_interpolate(field.at(y3, x0), field.at(y3, x1), field.at(y3, x2), field.at(y3, x3), x);

    cubic_interpolate(q0, q1, q2, q3, y)
}