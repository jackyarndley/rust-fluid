use crate::grid::Grid;

pub fn min<T: PartialOrd>(a: T, b: T) -> T {
    if a < b {
        return a
    }
    return b
}

pub fn max<T: PartialOrd>(a: T, b: T) -> T {
    if a > b {
        return a
    }
    return b
}

pub fn lerp(a: f64, b: f64, c: f64) -> f64 {
    a * (1.0 - c) + b * c
}

pub fn cerp(a: f64, b: f64, c: f64, d: f64, x: f64) -> f64 {
    let x_sq = x * x;
    let x_cu = x * x * x;

    let min_v = min(a, min(b, min(c, d)));
    let max_v = max(a, max(b, max(c, d)));

    let t = a * (0.0 - 0.5 * x + 1.0 * x_sq - 0.5 * x_cu) +
        b * (1.0 + 0.0 * x - 2.5 * x_sq + 1.5 * x_cu) +
        c * (0.0 + 0.5 * x + 2.0 * x_sq - 1.5 * x_cu) +
        d * (0.0 + 0.5 * x - 0.5 * x_sq + 0.5 * x_cu);

    min(max(t, min_v), max_v)
}

pub fn lerp_grid(grid: &Grid, x: f64, y: f64) -> f64 {
    let mut temp_x = min(max(x - grid.offset_x, 0.0), grid.columns as f64 - 1.001);
    let mut temp_y = min(max(y - grid.offset_y, 0.0), grid.rows as f64 - 1.001);

    let ix = temp_x as usize;
    let iy = temp_y as usize;

    temp_x -= ix as f64;
    temp_y -= iy as f64;

    let x00: f64 = grid.at(iy + 0, ix + 0);
    let x10: f64 = grid.at(iy + 0, ix + 1);
    let x01: f64 = grid.at(iy + 1, ix + 0);
    let x11: f64 = grid.at(iy + 1, ix + 1);

    lerp(lerp(x00, x10, temp_x), lerp(x01, x11, temp_x), temp_y)
}

pub fn cerp_grid(grid: &Grid, x: f64, y: f64) -> f64 {
    let mut temp_x = min(max(x - grid.offset_x, 0.0), grid.columns as f64 - 1.001);
    let mut temp_y = min(max(y - grid.offset_y, 0.0), grid.rows as f64 - 1.001);

    let ix = temp_x as usize;
    let iy = temp_y as usize;

    temp_x -= ix as f64;
    temp_y -= iy as f64;

    let x0 = max(ix as isize - 1, 0) as usize;
    let x1 = ix;
    let x2 = ix + 1;
    let x3 = min(ix + 2, grid.columns - 1);

    let y0 = max(iy as isize - 1, 0) as usize;
    let y1 = iy;
    let y2 = iy + 1;
    let y3 = min(iy + 2, grid.rows - 1);

    let q0: f64 = cerp(grid.at(y0, x0), grid.at(y0, x1), grid.at(y0, x2), grid.at(y0, x3), temp_x);
    let q1: f64 = cerp(grid.at(y1, x0), grid.at(y1, x1), grid.at(y1, x2), grid.at(y1, x3), temp_x);
    let q2: f64 = cerp(grid.at(y2, x0), grid.at(y2, x1), grid.at(y2, x2), grid.at(y2, x3), temp_x);
    let q3: f64 = cerp(grid.at(y3, x0), grid.at(y3, x1), grid.at(y3, x2), grid.at(y3, x3), temp_x);

    cerp(q0, q1, q2, q3, temp_y)
}

pub fn length(x: f64, y: f64) -> f64 {
    (x * x + y * y).sqrt()
}

pub fn cubic_pulse(x: f64) -> f64 {
    let y = min(x.abs(), 1.0);
    1.0 - y * y * (3.0 - 2.0 * y)
}