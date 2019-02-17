use crate::field::Field;

// Clamps a value between two bounds
pub fn clamp<T: PartialOrd>(value: T, lower: T, upper: T) -> T {
    if value < lower {
        lower
    } else if value > upper {
        upper
    } else {
        value
    }
}

pub fn linear_interpolate(a: f64, b: f64, c: f64) -> f64 {
    a * (1.0 - c) + b * c
}

pub fn empty(_: f64, _: f64, _: &Field) -> f64 {
    0.0
}

pub fn bilinear_interpolate(mut x: f64, mut y: f64, field: &Field) -> f64 {
    x -= field.offset_x;
    y -= field.offset_y;

    x = clamp(x, 0.0, field.columns as f64 - 1.001);
    y = clamp(y, 0.0, field.rows as f64 - 1.001);

    let p1_x = x.floor() as usize;
    let p1_y = y.floor() as usize;

    x -= p1_x as f64;
    y -= p1_y as f64;

    let p1 = field.at(p1_y, p1_x);
    let p2 = field.at(p1_y, p1_x + 1);
    let p3 = field.at(p1_y + 1, p1_x);
    let p4 = field.at(p1_y + 1, p1_x + 1);

    let l1 = linear_interpolate(p1, p2, x);
    let l2 = linear_interpolate(p3, p4, x);

    linear_interpolate(l1, l2, y)
}