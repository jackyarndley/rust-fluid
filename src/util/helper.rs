// Returns length of vector
pub fn length(x: f64, y: f64) -> f64 {
    (x * x + y * y).sqrt()
}

// Cubic pulse function returns value in range 0 - 1
pub fn cubic_pulse(mut x: f64) -> f64 {
    x = min(x.abs(), 1.0);
    1.0 - x * x * (3.0 - 2.0 * x)
}

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

pub fn min<T: PartialOrd>(a: T, b: T) -> T {
    if a < b {
        a
    } else {
        b
    }
}

pub fn max<T: PartialOrd>(a: T, b: T) -> T {
    if a > b {
        a
    } else {
        b
    }
}

pub fn nsgn(a: f64) -> f64 {
    if a < 0.0 {
        -1.0
    } else {
        1.0
    }
}

pub fn rotate(x: &mut f64, y: &mut f64, phi: f64) {
    let tmp_x = *x;
    let tmp_y = *y;

    *x = phi.cos() * tmp_x + phi.sin() * tmp_y;
    *y = -phi.sin() * tmp_x + phi.cos() * tmp_y;
}