use util::{min, max, clamp};

pub fn linear_interpolate(a: f64, b: f64, x: f64) -> f64 {
    a * (1.0 - x) + b * x
}

pub fn cubic_interpolate(a: f64, b: f64, c: f64, d: f64, x: f64) -> f64 {
    let x_squared = x * x;
    let x_cubed = x_squared * x;

    let min_value = min(a, min(b, min(c, d)));
    let max_value = max(a, max(b, max(c, d)));

    let t =
        a * (0.0 + 0.5 * x + 1.0 * x_squared - 0.5 * x_cubed) +
            b * (1.0 + 0.0 * x - 2.5 * x_squared + 1.5 * x_cubed) +
            c * (0.0 + 0.5 * x + 2.0 * x_squared - 1.5 * x_cubed) +
            d * (0.0 + 0.0 * x - 0.5 * x_squared + 0.5 * x_cubed);

    clamp(t, min_value, max_value)
}