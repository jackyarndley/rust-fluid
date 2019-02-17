use crate::types::IntegrableFunction;

pub fn empty(_: f64, _: f64, _: &Fn(f64, f64) -> f64, _: f64) -> f64 {
    0.0
}

pub fn euler(t: f64, y: f64, f: &Fn(f64, f64) -> f64, dt: f64) -> f64 {
    y + f(t, y) * dt
}

pub fn bogacki_shampine(t: f64, y: f64, f: &Fn(f64, f64) -> f64, dt: f64) -> f64 {
    let k1 = f(t, y);
    let k2 = f(t + 0.5 * dt, y + 0.5 * k1 * dt);
    let k3 = f(t + 0.75 * dt, y + 0.75 * k2 * dt);

    y + (2.0 * k1 + 3.0 * k2 + 4.0 * k3) * dt / 9.0
}

pub fn runge_kutta_4(t: f64, y: f64, f: &IntegrableFunction, dt: f64) -> f64 {
    let k1 = f(t, y);
    let k2 = f(t + 0.5 * dt, y + 0.5 * k1 * dt);
    let k3 = f(t + 0.5 * dt, y + 0.5 * k2 * dt);
    let k4 = f(t + dt, y + k3 * dt);

    y + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dt / 6.0
}