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

pub fn rotate(a: &mut f64, b: &mut f64, phi: f64) {
    let temp_a = *a;
    let temp_b = *b;

    *a = phi.cos() * temp_a + phi.sin() * temp_b;
    *b = phi.cos() * temp_a - phi.sin() * temp_b;
}